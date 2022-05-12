# -*- coding: utf-8 -*-

import requests
import pandas as pd
from tqdm import tqdm
from pubchempy import get_compounds


def normalize_chemicals():
    chemical_list = set(pd.read_csv(
        'data/normalized_data/graph_data.csv', sep='\t', usecols=['source']
    )['source'].to_list())

    chemical_data = []
    count = 0

    for chemical in tqdm(chemical_list):
        try:
            pubchem_compound = get_compounds(chemical, 'name')
            if len(pubchem_compound) > 0:
                count += 1
                pubchem_id = pubchem_compound[0].cid
                smiles = pubchem_compound[0].canonical_smiles
            else:
                pubchem_id = ''
                smiles = ''
        except IndexError:
            print(get_compounds(chemical, 'name'))
            raise ValueError

        chemical_data.append({
            'id': chemical,
            'pubchem id': pubchem_id,
            'smiles': smiles
        })

    print(f'Only {count} chemicals could be mapped to PubChem.')

    chemical_df = pd.DataFrame(chemical_data)
    chemical_df.to_csv('data/normalized_data/chemicals.csv', sep='\t', index=False)


def normalize_proteins():
    # TODO: Find another way to do this
    protein_list = set(pd.read_csv(
        'data/normalized_data/graph_data.csv', sep='\t', usecols=['target']
    )['target'].to_list())

    protein_data = []

    count = 0

    for num, protein in tqdm(enumerate(protein_list)):
        result = requests.post('http://grounding.indra.bio/ground', json={'text': protein}).json()

        if len(result) > 1:
            db = result[0]['term']['db'].lower()
            if db in ['hgnc', 'uniprot', 'fplx']:
                id = f'{db}:{result[0]["term"]["id"]}'
                count += 1
            else:
                id = ''
        else:
            id = ''

        uniqe_idx = str(num).zfill(4)

        protein_data.append({
            'id': f'PROT{uniqe_idx}',
            'name': protein,
            'identifier': id
        })

    print(f'{count} proteins were mapped to a database.')

    protein_df = pd.DataFrame(protein_data)
    protein_df.to_csv('data/normalized_data/protein.csv', sep='\t', index=False)


if __name__ == '__main__':
    normalize_proteins()