from preprocessor import resolver
import pandas as pd
from chembl_webresource_client.new_client import new_client
import requests
from bs4 import BeautifulSoup
from rdkit import DataStructs, Chem
from rdkit.Chem.Fingerprints import FingerprintMols
import numpy as np
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
from pubchempy import get_compounds


# Return a soup of the DrugBank page specified.
def soupify_page(db_cid):
    base_url = 'https://www.drugbank.ca/drugs/'
    abs_url = base_url + db_cid
    resp = requests.get(abs_url)
    return BeautifulSoup(resp.text, 'html5lib')


# Scrape and return the drug's SMILE to the user.
def get_db_smiles(db_cid):
    warm_soup = soupify_page(db_cid)
    text_to_find = 'SMILES'
    extr = warm_soup.find('dt', text=text_to_find)

    smile = None
    # Check if field exists for the specific drug.
    if extr:
        smile_acquired = extr.find_next_sibling('dd').text.strip()
        # If SMILE isn't available, return 'None'.
        if not smile_acquired == 'Not Available':
            smile = smile_acquired

    return smile


# Filter entity results by ChEMBL CID and return in DataFrame format.
def filter_results(entity, chembl_cid):
    # Overwrite any ChEMBL CIDs that are known to be bad.
    if chembl_cid == 'CHEMBL2029132':
        chembl_cid = 'CHEMBL2367706'

    results = getattr(new_client, entity)
    filtered_results = results.filter(molecule_chembl_id=chembl_cid)
    return pd.DataFrame(filtered_results)


# Return the Tanimoto similarity matrix given a DataFrame of SMILES.
def tanimoto_similarity(dataframe):
    smiles = list(dataframe.smile)

    # Make a list of RDKit fingerprints first.
    fps = []
    for smile in smiles:
        fps_found = None
        # Skip over any loop steps that involve SMILES of type None.
        if smile:
            mol = Chem.MolFromSmiles(smile)
            # And, again, skip any molecules that couldn't be found.
            if mol:
                fps_found = FingerprintMols.FingerprintMol(mol)
        fps.append(fps_found)

    # Remove any drugs without a fingerprint?
    fps = [each for each in fps if each is not None]

    num_of_drugs = len(fps)
    print(f'{len(smiles)-num_of_drugs} drugs removed')

    # Matrix is symmetric, thus only computing the upper-triangular part should suffice.
    upper_triangular = np.zeros(int((num_of_drugs**2-num_of_drugs)/2))
    for _idx, pair in enumerate(combinations(enumerate(fps), 2)):
        if pair[0][1] and pair[1][1]:
            dist = DataStructs.FingerprintSimilarity(pair[0][1], pair[1][1],
                                                     metric=DataStructs.TanimotoSimilarity)
        # If any of the fingerprints doesn't exist, set coefficient to zero.
        else:
            dist = 0
        # Append to the flattened upper triangular matrix.
        upper_triangular[_idx] = dist

    tanimoto_matrix = np.zeros((num_of_drugs, num_of_drugs))

    tri_upper_indices = np.triu_indices(num_of_drugs, k=1)
    tanimoto_matrix[tri_upper_indices] = upper_triangular
    tanimoto_matrix = tanimoto_matrix + tanimoto_matrix.T
    np.fill_diagonal(tanimoto_matrix, 1)

    return tanimoto_matrix


if __name__ == '__main__':
    [modified_sheets, _] = resolver()
    for sheet in modified_sheets:
        acquired_data = {'smile': []}
        for idx, drug in enumerate(sheet.itertuples()):
            # Set DrugBank as preferred source...
            drug_smile = get_db_smiles(drug.DrugBank_CID)

            # ...then default to ChEMBL when no SMILES are returned...
            if not drug_smile:
                drug_info = filter_results('molecule', drug.ChEMBL_CID)
                molecule_structures = drug_info.molecule_structures[0]
                if molecule_structures:
                    drug_smile = molecule_structures['canonical_smiles']
                # ...and finally prod PubChem's API as last resort.
                else:
                    drug_name = drug_info.pref_name[0]
                    pubchem_info = get_compounds(drug_name, 'name')
                    if pubchem_info:
                        drug_smile = pubchem_info[0].canonical_smiles
                    else:
                        # If no SMILE's found, return None.
                        drug_smile = None

            acquired_data['smile'].append(drug_smile)

        # Append all new data to the original DataFrame!
        column = pd.DataFrame.from_dict(acquired_data)
        _sheet = pd.concat([sheet, column], axis=1)

    similarity_matrix = tanimoto_similarity(_sheet)
    # df = pd.DataFrame(data=similarity_matrix,
    #                   index=_sheet.DrugBank_CID,
    #                   columns=_sheet.ChEMBL_CID)
    df = pd.DataFrame(data=similarity_matrix)

    # Display and save a cluster-map of the results.
    sns.clustermap(df, yticklabels=False, xticklabels=False)
    plt.savefig('tanimoto_clustermap.png', dpi=300,)
    plt.show()

    # Save smiles and tanimoto matrix to spreadsheet.
    _sheet.to_excel('output_smiles.xlsx')
    df.to_excel('output_tanimoto_similarity.xlsx')
