import pandas as pd
import pathlib
import os
from typing import Union
import re
import requests
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support import expected_conditions as ec
from selenium.webdriver.support.ui import WebDriverWait


# Find and import every spreadsheet found in the script path.
def _importer():
    script_path = pathlib.Path(__file__).parent.absolute()

    # PyCharm BUG: pathlib.Path objects not considered 'PathLike'
    fn: Union[pathlib.Path, os.PathLike]
    # Retrieve all Excel files, but skip any starting with 'output' or 'temp'.
    sheet_paths = [fn for fn in script_path.glob('*.xlsx')
                   if not os.path.basename(fn).startswith(('output', 'temp'))]

    sheet_list = []
    for path in sheet_paths:
        sheet = pd.read_excel(path)
        sheet_list.append(sheet)
    num_of_sheets = len(sheet_list)

    return sheet_list, num_of_sheets


# Use UniChem's REST API to get the associated DB/ChEMBL CID.
def _request_cid(source_cid):
    base_url = 'https://www.ebi.ac.uk/unichem/rest'
    # Check if source CID comes from DB or ChEMBL.
    if re.match('DB', source_cid):
        ext = '/2/1'
    else:
        ext = '/1/2'
    abs_url = base_url + '/src_compound_id/' + source_cid + ext
    resp = requests.get(abs_url)
    # Only proceed if an 'OK 200' response was given.
    if resp:
        resp_json = resp.json()
        # Make sure the answer is not an empty list.
        if resp_json:
            return resp_json[0]['src_compound_id']
        else:
            return None
    else:
        return None


# Simple wrapper for ChEMBL's REST API.
def prod_chembl_api(cid, query):
    # Request that the response is given in JSON format.
    headers = {
        'Accept': 'application/json',
        'Content-Type': 'application/json; charset=UTF-8'
    }
    chembl_url = 'https://www.ebi.ac.uk/chembl/api/data/molecule/search?q='
    resp = requests.get(chembl_url + cid, headers=headers)
    resp_json = resp.json()
    if resp_json:
        return resp_json['molecules'][0][query]
    else:
        return None


def _fget_chembl_cid(db_cid):
    base_url = 'https://www.drugbank.ca/drugs/'
    abs_url = base_url + db_cid
    resp = requests.get(abs_url)
    soup = BeautifulSoup(resp.text, 'html5lib')
    extr_dt = soup.find('dt', string='Name')
    # This might help with garbage collection...
    soup.decompose()
    extr = extr_dt.find_next_sibling('dd').text
    return prod_chembl_api(extr, 'molecule_chembl_id')


# Return ChEMBL CID when given a DrugBank CID as input.
def _get_chembl_cid(db_cid):
    base_url = 'https://www.drugbank.ca/drugs/'
    abs_url = base_url + db_cid
    resp = requests.get(abs_url)
    # Scrape the page for the ChEMBL link.
    soup = BeautifulSoup(resp.text, 'html5lib')
    extr = soup.find('dt', string='ChEMBL')
    soup.decompose()
    # Check if scraping was successful...
    if not extr:
        # If not, resort to UniChem.
        chembl_cid = _request_cid(db_cid)
        if not chembl_cid:
            # Or search and scrape!
            return _fget_chembl_cid(db_cid)
        else:
            return chembl_cid
    # If so, extract and return ChEMBL CID.
    else:
        extr = extr.find_next_sibling('dd')
        return extr.text


# Make it a class for easier Selenium webdriver garbage management.
class _DbCID:
    def __init__(self):
        # Initialise a headless chrome webdriver in selenium.
        script_path = pathlib.Path(__file__).parent.absolute()
        driver_path = '/chromedriver_win32/chromedriver.exe'
        abs_path = str(script_path) + driver_path
        chrome_options = Options()
        chrome_options.add_argument("--headless")
        self.driver = webdriver.Chrome(
            executable_path=abs_path, options=chrome_options)

    def __del__(self):
        self.driver.quit()

    def _fget(self, chembl_cid):
        name = prod_chembl_api(chembl_cid, 'pref_name')
        self.driver.get('https://www.drugbank.ca')
        # Search for the drug in question and click 'Enter'.
        search_locator = "(//input[@id='query'])"
        self.driver.find_element_by_xpath(search_locator).send_keys(name + '\ue007')
        current_url = self.driver.current_url
        db_cid = current_url.split('/')[-1]
        # Check if extracted string is indeed a DB CID...
        if re.match('DB', db_cid):
            answer = db_cid
        # ...otherwise we've landed in a search results page!
        else:
            # Click the first result, wait for page to load...
            # TODO: Match alternative and trade names before accepting the first result as correct.
            first_result_locator = "(//div[contains(@class,'unearth-drug-search-results')]//div//div//h2)"
            self.driver.find_element_by_xpath(first_result_locator).click()
            WebDriverWait(self.driver, 15).until(ec.url_changes(current_url))
            # ...and extract DB CID!
            db_cid = self.driver.current_url.split('/')[-1]
            if re.match('DB', db_cid):
                answer = db_cid
            else:
                answer = None
        # There's only one working webdriver tab open, no need to close it.
        # self.driver.close()
        return answer

    # Return ChEMBL CID when given a DrugBank CID as input.
    def get(self, chembl_cid):
        # Query UniChem first.
        db_cid = _request_cid(chembl_cid)
        if not db_cid:
            # And if it fails, search and scrape.
            return self._fget(chembl_cid)
        else:
            return db_cid


# Map and append ChEMBL to DrugBank CIDs and vice versa.
def resolver():
    # This object's responsible of the Chrome webdriver.
    db_cid = _DbCID()

    imported_sheets = _importer()
    modified_sheets = []
    num_of_drugs_per_sheet = []
    for sheet in imported_sheets[0]:
        # Comprehensive list of CIDs: [0]:DB, [1]:ChEMBL
        list_of_cids = [[] for _ in range(2)]
        drug_series = sheet.iloc[:, 1]

        # Only return 'true' for DrugBank CIDs.
        is_db_cid_series = drug_series.str.match('DB')
        for idx, is_db in enumerate(is_db_cid_series):
            drug_cid = drug_series[idx]
            if is_db:
                retrieved_cid = _get_chembl_cid(drug_cid)
                list_of_cids[0].append(drug_cid)
                list_of_cids[1].append(retrieved_cid)
            else:
                retrieved_cid = db_cid.get(drug_cid)
                list_of_cids[0].append(retrieved_cid)
                list_of_cids[1].append(drug_cid)
            print(str(idx) + ': ' + drug_cid + ' > ' + retrieved_cid)

        # Append CIDs to the DataFrame in a structured manner.
        sheet['DrugBank_CID'], sheet['ChEMBL_CID'] = [list_of_cids[0],
                                                      list_of_cids[1]]
        modified_sheets.append(sheet)
        num_of_drugs_per_sheet.append(len(sheet.index))
    # Destroying it, also destroys the Chrome webdriver!
    del db_cid
    # Return a list of modified sheets and their length.
    return modified_sheets, num_of_drugs_per_sheet
