# -*- coding: utf-8 -*-

import pandas as pd
from tqdm import tqdm


def preprocess_file():

    df = pd.read_csv(
        'data/raw/NMR-Protein-Screening-hits_Protein_Sequence.txt',
        sep='\t',
        dtype=str,
        encoding='unicode_escape',
        skiprows=[1,2],
        engine='python'
    )

    graph_data = []

    for prot_name in tqdm(df.columns):
        for idx in range(0, df.shape[0]):
            chem = df.loc[idx, prot_name]

            if pd.notna(chem):
                graph_data.append({
                    'source': chem,
                    'relation': 'binds',
                    'target': prot_name
                })

    graph_df = pd.DataFrame(graph_data)
    graph_df.to_csv('data/normalized_data/graph_data.csv', sep=',', index=False)


if __name__ == '__main__':
    preprocess_file()