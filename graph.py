# -*- coding: utf-8 -*-

"""Create BEL graph"""

import os
import pandas as pd
from tqdm import tqdm
from pubchempy import Compound

import pybel
from pybel import BELGraph
from pybel.dsl import Abundance, Protein, Entity

RELATION_MAPPER = {
    'binds': BELGraph.add_binds
}


def create_chemical_nodes(node_dict: dict):
    """Method to create chemical nodes"""
    df = pd.read_csv('data/normalized_data/chemicals.csv', sep='\t')

    for row in tqdm(df.values, desc='Creating chemical nodes'):
        (
            id,
            pubchem_id,
            smiles
        ) = row

        if pd.notna(pubchem_id):
            pubchem_id = int(pubchem_id)  # convert flot to int values

            chemical_node = Abundance(
                namespace='pubchem.compound',
                identifier=str(pubchem_id),
                name=Compound.from_cid(pubchem_id).synonyms[0],
                xrefs=[Entity(
                    namespace='enamine',
                    identifier=id,
                ), Entity(
                    namespace='smiles',
                    identifier=smiles,
                )]
            )
        else:
            chemical_node = Abundance(
                namespace='enamine',
                identifier=id,
            )

        node_dict[id] = chemical_node

    return node_dict


def create_protein_nodes(node_dict: dict):
    """Method to create protein nodes"""
    df = pd.read_csv('data/normalized_data/protein.csv', sep='\t')

    for row in tqdm(df.values, desc='Creating protein nodes'):
        (
            id,
            name,
            uniprot
        ) = row

        if pd.notna(uniprot):
            protein_node = Protein(
                namespace='uniprot',
                identifier=uniprot,
                name=name,
                xrefs=[Entity(
                    namespace='prot',
                    identifier=id,
                )]
            )
        else:
            protein_node = Abundance(
                namespace='prot',
                identifier=id,
            )

        node_dict[name] = protein_node

    return node_dict

def create_chembl_nodes(node_dict: dict):
    """Method to create protein nodes"""
    df = pd.read_excel('data/normalized_data/protein.csv', sep='\t')



def save_graph(bel_graph: BELGraph):
    """Save the BEL graph"""
    os.makedirs('data/graph', exist_ok=True)
    return pybel.dump(bel_graph, path='data/graph/covid_nmr.bel.nodelink.json')


def create_graph():
    """Method to create BEL graph from NMR data"""
    nodes = {}

    nodes = create_chemical_nodes(node_dict=nodes)
    nodes = create_protein_nodes(node_dict=nodes)

    graph_df = pd.read_csv('data/normalized_data/graph_data.csv', sep='\t')

    graph = BELGraph(name='COVID-NMR')

    for row in tqdm(graph_df.values, desc='Exporting data to BEL'):
        (
            source_id,
            rel,
            target_id
        ) = row
        source_node = nodes[source_id]
        target_node = nodes[target_id]

        RELATION_MAPPER[rel](
            graph,
            source_node, target_node,
            evidence='From COVID_NMR data',
            citation=('database', 'covid-nmr')
        )

    graph.summarize()

    save_graph(graph)


if __name__ == '__main__':
    create_graph()


