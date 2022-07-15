import pybel
import pandas as pd
from pybel.dsl import Protein
from pybel.dsl import Abundance
from pybel.dsl import Pathology
from pybel.dsl import BiologicalProcess
from pybel.dsl import Population
from pybel.dsl import Gene
from pybel.dsl import MicroRna
from pybel.dsl import Rna
from pybel.dsl import Fragment
import chembl_webresource_client
import openpyxl
import networkx as nx
from pybel.io.jupyter import to_jupyter
import matplotlib.pyplot as plt
import chembl_webresource_client
from chembl_webresource_client.new_client import new_client
import pubchempy
import pickle
import re
from urllib.parse import urlparse
import urllib
import requests
import io

#function to retrieve mechanisms from ChEMBL
#nmrGraph = pybel.BELGraph(name='NMRgraph')
##nmrData = pd.read_csv('NMR_data_new.csv')
##nmr2chembl = pd.read_excel('Protein-Screening-Summary_CHEMBL30_Actives_teams_max08tan.xlsx')
#itmp_chem = pd.read_csv('C:\\Users\\reagon.karki\\Documents\\ITMP\\EOSC Future\\ITMP ChEMBL data\\CHEMBL4495564.csv',sep=';',usecols=['ChEMBL ID'])
##nmr2pchem = pd.read_excel('PubChem_Actives_Literature_NCATS.xlsx',sheet_name='Similar_065')

#Function to retrieve mechanism of actions and target proteins from ChEMBL
#Returns a dictionary
def RetMech(chemblIds):
    getMech = new_client.mechanism
    mechList = []
    for i in range(len(chemblIds)):
        mechs = getMech.filter(molecule_chembl_id=chemblIds[i]).only(['mechanism_of_action','target_chembl_id'])
        #mechs = getMech.filter(molecule_chembl_id=chemblIds[i])
        print(mechs)
        mechList.append(list(mechs))
    named_mechList = dict(zip(chemblIds,mechList))
    named_mechList = {k: v for k, v in named_mechList.items() if v}
    return(named_mechList)

#Function to retrieve associated diseases
#Returns a dictionary
def RetDrugInd(chemblIDs):
    getDrugInd = new_client.drug_indication
    drugIndList = []
    for i in range(len(chemblIDs)):
        drugInd = getDrugInd.filter(molecule_chembl_id=chemblIDs[i]).only('mesh_heading')
        #drugInd = getDrugInd.filter(molecule_chembl_id=chemblIDs[i])
        print(drugInd)
        drugIndList.append(list(drugInd))
    named_drugIndList = dict(zip(chemblIDs,drugIndList))
    named_drugIndList = {k: v for k, v in named_drugIndList.items() if v}
    return(named_drugIndList)

#Function to retrieve associated assays
#Returns a dictionary
def RetAct(chemblIds,j):
    GetAct = new_client.activity
    ActList = []
    #for i in range(len(chemblIds)):
    #print(chemblIds[0])
    for i in range(len(chemblIds)):
        filtered_list=['assay_chembl_id','assay_type','pchembl_value','target_chembl_id','target_organism','bao_label']
        acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False,assay_type_iregex='(B|F)',target_organism='Homo sapiens').only(filtered_list)
        #acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False,target_organism='Homo sapiens').only(filtered_list)
        #acts = GetAct.filter(molecule_chembl_id=chemblIds[i],pchembl_value__isnull=False)
        j=j+1

        #acts = [d for d in acts if d.get('target_organism') == 'Homo sapiens']
        #acts = [d for d in acts if d.get('bao_label') == 'single protein format']
        #acts = [d for d in acts if d.get('type') in ['Ki', 'IC50']]
        acts = [d for d in acts if float(d.get('pchembl_value')) >= 6]
        print(j)
        #print(len(acts))
        acts = acts[:5]
        print(acts)
        ActList.append(list(acts))
    #print(ActList)
    named_ActList = dict(zip(chemblIds,ActList))
    named_ActList = {k: v for k, v in named_ActList.items() if v}
    return(named_ActList)

def chembl2uniprot(chemblIDs, count):
    getTarget = new_client.target
    chem2Gene2path = []
    for i in range(len(chemblIDs)):
        print(count)
        count = count + 1
        chem2path = []
        chem = getTarget.filter(chembl_id=chemblIDs[i]).only('target_components')

        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]
        print(getGene)

        chem2path = [item for item in chem[0]['target_components'][0]['target_component_xrefs'] if
                     item["xref_src_db"] == "Reactome"]

        chem2path.append(getGene)
        chem2Gene2path.append(list(chem2path))

    named_chem2Gene2path = dict(zip(chemblIDs, chem2Gene2path))
    named_chem2Gene2path = {k: v for k, v in named_chem2Gene2path.items() if v}
    return (named_chem2Gene2path)

def chembl2gene2path(chem2geneList,ActList):
    for item in chem2geneList:
        #print(item)
        sizeOfitem = len(chem2geneList[item])
        #print(chem2Gene[item])
        #print(chem2geneList[item][sizeOfitem-1]['component_synonym'])
        gene = chem2geneList[item][sizeOfitem-1]['component_synonym']
        for jtem in ActList:
            #print(jtem,'b')
            #print(len(ActList_copy[jtem]))
            for i in range(len(ActList[jtem])):
                if item == ActList.get(jtem)[i]['target_chembl_id']:
                #print(jtem)
                    newkey = {'Protein': gene}
                    ActList[jtem][i].update(newkey)
                #print(ActList_copy.get(jtem)[i]['target_chembl_id'])

    return(ActList)


def filter_graph(mainGraph, vprotList):
    nsp_list = []
    for u, v, data in mainGraph.edges(data=True):
        if u.name in vprotList or v.name in vprotList:
            nsp_list.append(u)
            nsp_list.append(v)

    for u, v, data in mainGraph.edges(data=True):
        if 'Similarity' not in data:
            continue
        if data['Similarity'] >= .65 and (u in nsp_list or v in nsp_list):
            # yield u, v, k, data
            nsp_list.append(u)
            nsp_list.append(v)

    # print(nsp_list)

    nsp_graph = mainGraph.subgraph(nsp_list)
    # return(to_jupyter(nsp_graph))
    final_nsp = [node for node in nsp_graph.nodes() if isinstance(node, pybel.dsl.Protein) and node.name in vprotList]
    # print(final_nsp)
    for u, v, data in nsp_graph.edges(data=True):
        if u.namespace in ['Enamine', 'ChEMBL', 'CID'] and v.namespace in ['Enamine', 'ChEMBL', 'CID']:
            final_nsp.append(u)
            final_nsp.append(v)

    # return(nsp_graph.subgraph(final_nsp))
    nsp_graph = nsp_graph.subgraph(final_nsp)
    # chembl_act = []
    nsp_nodes = [node for node in nsp_graph.nodes()]
    for u in nsp_graph.nodes():
        if u.namespace == 'ChEMBL':
            # test = nmrGraph.neighbors(u)
            chem = [n for n in mainGraph.neighbors(u)]
            nsp_nodes.extend(chem)

    return (mainGraph.subgraph(nsp_nodes))


#Functions for creating graph
def chem2moa_rel(named_mechList,itmpGraph):
    for i in named_mechList:
    #print(i)
    #break
    #print(named_mechList[i])
    #print(len(named_mechList[i]))
    #break
        for j in range(len(named_mechList[i])):
            #print(named_mechList[i][j]['mechanism_of_action'])
            #print(named_mechList[i][j]['target_chembl_id'])
            #print(i)
            #break
            itmpGraph.add_association(Abundance(namespace='ChEMBL',name=i),BiologicalProcess(namespace='MOA',name=named_mechList[i][j]['mechanism_of_action']),citation='ChEMBL database',evidence='ChEMBL query')
            if not named_mechList[i][j]['target_chembl_id'] == None:
                itmpGraph.add_association(Abundance(namespace='ChEMBL',name=i),Protein(namespace='HP',name=named_mechList[i][j]['Protein']),citation='ChEMBL database',evidence='ChEMBL query')
    return(itmpGraph)

def chem2dis_rel(named_drugIndList,itmpGraph):
    for i in named_drugIndList:
    #print(i)
    #break
    #print(named_drugIndList[i])
    #print(len(named_drugIndList[i]))
    #break
        for j in range(len(named_drugIndList[i])):
            #print(named_drugIndList[i][j]['mesh_heading'])
            #print(i)
            #break
            itmpGraph.add_association(Abundance(namespace='ChEMBL',name=i),Pathology(namespace='Disease',name=named_drugIndList[i][j]['mesh_heading']),citation='ChEMBL database',evidence='ChEMBL query')
    return(itmpGraph)

def chem2act_rel(named_ActList,itmpGraph):
    for i in named_ActList:
        for j in range(len(named_ActList[i])):
            # print(named_mechList[i][j]['mechanism_of_action'])
            # print(named_mechList[i][j]['target_chembl_id'])
            # print(i)
            # break
            # nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),BiologicalProcess(namespace='MOA',name=mechList[i][j]['mechanism_of_action']),citation='ChEMBL database',evidence='ChEMBL query')
            if not named_ActList[i][j]['target_chembl_id'] == None:
                itmpGraph.add_association(Abundance(namespace='ChEMBLAssay',name=named_ActList[i][j]['assay_chembl_id']),
                                         Protein(namespace='HP', name=named_ActList[i][j]['ProteinName']),
                                         citation='ChEMBL database', evidence='ChEMBL query')

            itmpGraph.add_association(Abundance(namespace='ChEMBL', name=i),Abundance(namespace='ChEMBLAssay',name=named_ActList[i][j]['assay_chembl_id']),
                                      citation='ChEMBL database', evidence='ChEMBL query',assayType=named_ActList[i][j]['assay_type'],
                                      pChEMBL=named_ActList[i][j]['pchembl_value'])

    return(itmpGraph)

def chem2act_rel_2(named_ActList, itmpGraph):
    for i in named_ActList:
        # print(i)
        for j in range(len(named_ActList[i])):
            # print(named_mechList[i][j]['mechanism_of_action'])
            # print(named_mechList[i][j]['target_chembl_id'])
            # print(i)
            # break
            # nmrGraph.add_association(Pathology(namespace='ChEMBL',name=i),BiologicalProcess(namespace='MOA',name=mechList[i][j]['mechanism_of_action']),citation='ChEMBL database',evidence='ChEMBL query')
            if not named_ActList[i][j]['target_chembl_id'] == None:
                if 'Protein' in named_ActList[i][j]:
                    itmpGraph.add_association(
                        Abundance(namespace='ChEMBLAssay', name=named_ActList[i][j]['assay_chembl_id']),
                        Protein(namespace='HP', name=named_ActList[i][j]['Protein']),
                        citation='ChEMBL database', evidence='ChEMBL query')
                else:
                    itmpGraph.add_association(
                        Abundance(namespace='ChEMBLAssay', name=named_ActList[i][j]['assay_chembl_id']),
                        Protein(namespace='HP', name=named_ActList[i][j]['target_chembl_id']),
                        citation='ChEMBL database', evidence='ChEMBL query')

            itmpGraph.add_association(Abundance(namespace='ChEMBL', name=i),
                                      Abundance(namespace='ChEMBLAssay', name=named_ActList[i][j]['assay_chembl_id']),
                                      citation='ChEMBL database', evidence='ChEMBL query',
                                      assayType=named_ActList[i][j]['assay_type'],
                                      pChEMBL=named_ActList[i][j]['pchembl_value'])

    return (itmpGraph)

# def chem2gene2path_rel(named_chem2geneList,itmpGraph):
#     for item in named_chem2geneList:
#         itemLen = len(named_chem2geneList[item])-2
#         for j in range(itemLen):
#             #print(named_chem2geneList)
#             itmpGraph.add_association(MicroRna(namespace='HP', name=named_chem2geneList[item][itemLen]['component_synonym']),
#                                       bp(namespace='Pathway',name=named_chem2geneList[item][j]['xref_name']),
#                                       citation='ChEMBL database', evidence='ChEMBL query',
#                                       Reactome=named_chem2geneList[item][j]['xref_id'])
#
#     return(itmpGraph)


def chembl2uniprot(chemblIDs, count):
    getTarget = new_client.target
    chem2Gene2path = []
    chemHasNoPath = []
    chemNotprotein = []
    for i in range(len(chemblIDs)):
        # print(count)
        count = count + 1
        chem2path = []
        chem = getTarget.filter(chembl_id=chemblIDs[i]).only('target_components')
        #print(chem)
        # break
        try:
            uprot_id = chem[0]['target_components'][0]['accession']
        except IndexError:
            # print(chemblIDs[i])
            chemHasNoPath.append(chemblIDs[i])
            continue

        if chem[0]['target_components'][0]['accession'] == None:
            chemHasNoPath.append(chemblIDs[i])
            # continue

    chemblIDs_clean = [item for item in chemblIDs if item not in chemHasNoPath]
    print('old', len(chemblIDs_clean))
    for i in range(len(chemblIDs_clean)):

        chem = getTarget.filter(chembl_id=chemblIDs_clean[i]).only('target_components')
        # print(chem)
        # break
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        # rint(getGene)
        # break
        try:
            getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]
        except IndexError:
            chemNotprotein.append(chemblIDs_clean[i])
            continue
    chemblIDs_clean = [item for item in chemblIDs_clean if item not in chemNotprotein]
    print('newLen', len(chemblIDs_clean))
    print(len(chemNotprotein))
    # break
    # print(getGene)
    # break
    for i in range(len(chemblIDs_clean)):
        chem = getTarget.filter(chembl_id=chemblIDs_clean[i]).only('target_components')
        # print(chem)
        # break
        uprot_id = chem[0]['target_components'][0]['accession']
        
        getGene = chem[0]['target_components'][0]['target_component_synonyms']
        #print(getGene)
        getGene = [item for item in getGene if item["syn_type"] == "GENE_SYMBOL"][0]
        #print(getGene)

        chem2path = [item for item in chem[0]['target_components'][0]['target_component_xrefs'] if
                     item["xref_src_db"] == "Reactome"]
        
        uprot = {'accession':uprot_id}
        #print(uprot)
        #print(getGene)
        chem2path.append(uprot)
        chem2path.append(getGene)
        #print(chem2path)
        # break
        chem2Gene2path.append(chem2path)
        
        # print(chem2Gene2path)
        # break

        # print(chem2Gene2path)

    # print(chemHasNoPath)
    named_chem2Gene2path = dict(zip(chemblIDs_clean, chem2Gene2path))
    named_chem2Gene2path = {k: v for k, v in named_chem2Gene2path.items() if v}
    return (named_chem2Gene2path)

# Functions for creating graph
def chem2moa_rel_2(named_mechList, itmpGraph):
    for i in named_mechList:
        # print(i)
        # break
        # print(named_mechList[i])
        # print(len(named_mechList[i]))
        # break
        for j in range(len(named_mechList[i])):
            # print(named_mechList[i][j]['mechanism_of_action'])
            # print(named_mechList[i][j]['Protein'])
            # print(i)
            # break
            itmpGraph.add_association(Abundance(namespace='ChEMBL', name=i), BiologicalProcess(namespace='MOA', name=
            named_mechList[i][j]['mechanism_of_action']), citation='ChEMBL database', evidence='ChEMBL query')
            if not named_mechList[i][j]['target_chembl_id'] == None:
                # print(named_mechList[i])
                if 'Protein' in named_mechList[i][j]:
                    # print('yes')
                    itmpGraph.add_association(Abundance(namespace='ChEMBL', name=i),
                                              Protein(namespace='HP', name=named_mechList[i][j]['Protein']),
                                              citation='ChEMBL database', evidence='ChEMBL query')
                else:
                    itmpGraph.add_association(Abundance(namespace='ChEMBL', name=i),
                                              Protein(namespace='HP', name=named_mechList[i][j]['target_chembl_id']),
                                              citation='ChEMBL database', evidence='ChEMBL query')

    return (itmpGraph)

def Ret_chembl_protein(sourceList):
    protein_List = []
    for item in sourceList:
        for j in range(len(sourceList[item])):
            protein_List.append(sourceList[item][j]['target_chembl_id'])

    protein_List = set(protein_List)
    protein_List = list(filter(None, protein_List))
    return (protein_List)

# def ExtractFromUniProt(uniprot_id):
    # from bioservices import UniProt
    # Uniprot_Dict = []
    # # Make a link to the UniProt webservice
    # service = UniProt()

    # for id in uniprot_id:

        # # create URL for each uniprot id
        # url = 'https://www.uniprot.org/uniprot/' + id + '.txt'
        # print(url)

        # # Retrieve data for id in text format
        # ret_uprot = urllib.request.urlopen(url)

        # print(id)
        # id_copy = id
        # i = 0
        # j = 0
        # id = {}
        # id['Disease'] = {}
        # id['Reactome'] = {}
        # id['Function'] = {}
        # id['BioProcess'] = {}
        # # print(id)

        # # parse each line looking for info about disease, pathway, funcn, bp and so on
        # for line in ret_uprot:

            # line = line.decode('utf-8')

            # # parse lines with disease and extract disease names and omim ids
            # if '-!- DISEASE:' in line:
                # if ('[MIM:' in line):
                    # dis = line.split(':')
                    # # dis returns list of splitted text, [1] = name of dis, [2] = OMIM ID, extra chars need cleaning
                    # # print(dis[1][1:-5])
                    # # print(dis[2][:-1])
                    # id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})

            # # extract reactome ids and names
            # if 'Reactome;' in line:
                # ract = line.split(';')
                # # ract returns list with reactome id and name, needs cleaning
                # id['Reactome'].update({ract[2][1:-2]: ract[1][1:]})
                # # print(ract[1][1:])
                # # print(ract[2][1:-2])

            # # look for functions
            # if ' F:' in line:
                # if j < 5:
                    # # take only first 5 entries for now
                    # # print(j)
                    # fn = line.split(';')
                    # # fn returns list with GO ids and names
                    # id['Function'].update({fn[2][3:]: fn[1][1:]})
                    # # print(fn[1][1:])
                    # # print(fn[2][3:])
                    # j += 1

            # # look for biological processes
            # if ' P:' in line:
                # if i < 5:
                    # # take only first 5 entries for now
                    # # print(i)
                    # bp = line.split(';')
                    # # bp returns list with GO ids and names
                    # id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
                    # # print(bp[1][1:])
                    # # print(bp[2][3:])
                    # i += 1

        # # fetch info about gene, len, seq, mass from api directly
        # # Make a query string
        # query = "accession:" + str(id_copy)

        # # Define a list of columns we want to retrive
        # # columnlist = "id,entry name,length,mass,go(biological process),go(molecular function), pathway,feature(TOPOLOGICAL DOMAIN),comment(DISEASE),pathway"
        # columnlist = "id,entry name,genes(PREFERRED),length,mass,sequence"

        # # Run the remote search
        # result = service.search(query, frmt="tab", columns=columnlist)

        # df_result = pd.read_table(io.StringIO(result))

        # geneName = {'Gene': df_result['Gene names  (primary )'][0]}
        # seqLen = {'Length': df_result['Length'][0]}
        # # remove comma from mass values
        # mass = {'Mass': int(df_result['Mass'][0].replace(",", ""))}
        # seq = {'Sequence': df_result['Sequence'][0]}

        # id.update(geneName)
        # id.update(seqLen)
        # id.update(mass)
        # id.update(seq)

        # Uniprot_Dict.append(id)

    # Uniprot_Dict = dict(zip(uniprot_id, Uniprot_Dict))

    # return(Uniprot_Dict)
    
    
#function to convert pubchem cids to chembl ids

def cid2chembl(cidList):
    import pubchempy as pcp
    cid2chembl_list = []
    for id in cidList:
        #GetChembl = {}
        c=pcp.Compound.from_cid(id)
        syn = c.synonyms
        for s in syn:
            if s.startswith('CHEMBL'):
                #GetChembl['CHEMBL'] = s
                cid2chembl_list.append(s)
                print('Pubchem:',id,'Converted to:',s)

        #cid2chembl_list.append(GetChembl)
    return(cid2chembl_list) 
    
# def uniprot_rel(named_uprotList,itmpGraph):
    
    # for item in named_uprotList:
        # #print(named_uprotList[item]['Function'].keys())
        # fun=list(named_uprotList[item]['Function'].keys())
        # bp = list(named_uprotList[item]['BioProcess'].keys())
        # for f in fun:
            # if str(named_uprotList[item]['Gene']) != 'nan':
                # itmpGraph.add_association(Protein(namespace='MP',name=named_uprotList[item]['Gene']),MicroRna(namespace='Function',name=f),
                                         # citation='UniProt database',evidence='UniProt query')
            # else:
                # itmpGraph.add_association(Protein(namespace='MP',name=item),MicroRna(namespace='Function',name=f),
                                         # citation='UniProt database',evidence='UniProt query')
                
        
        # for b in bp:
            # if str(named_uprotList[item]['Gene']) != 'nan':
                # itmpGraph.add_association(Protein(namespace='MP',name=named_uprotList[item]['Gene']),BiologicalProcess(namespace='BioPro',name=b),
                                         # citation='UniProt database',evidence='UniProt query')
            # else:
                # itmpGraph.add_association(Protein(namespace='MP',name=item),MicroRna(namespace='Function',name=b),
                                         # citation='UniProt database',evidence='UniProt query')
        
    # return(itmpGraph)
    
# def ExtractFromUniProt(uniprot_id):
#     from bioservices import UniProt
#     Uniprot_Dict = []
#     # Make a link to the UniProt webservice
#     service = UniProt()
#
#     for id in uniprot_id:
#
#         # create URL for each uniprot id
#         url = 'https://www.uniprot.org/uniprot/' + id + '.txt'
#         print(url)
#
#         #Retrieve data for id in text format if found in uniprot
#
#         try:
#             ret_uprot = urllib.request.urlopen(url)
#         except urllib.request.HTTPError:
#             #uniprot_id.remove(id)
#             continue
#
#
# #         with urllib.request.urlopen(url) as response:
# #            ret_uprot = response.read()
#
#         print(id)
#         id_copy = id
#         i = 0
#         j = 0
#         k = 0
#         id = {}
#         id['Disease'] = {}
#         id['Reactome'] = {}
#         id['Function'] = {}
#         id['BioProcess'] = {}
#         id['Gene'] = {}
#         #id['Gene'] = {}
#         # print(id)
#
#         # parse each line looking for info about disease, pathway, funcn, bp and so on
#         for line in ret_uprot:
#
#             line = line.decode('utf-8')
#
#             # parse lines with disease and extract disease names and omim ids
#             if '-!- DISEASE:' in line:
#                 if ('[MIM:' in line):
#                     dis = line.split(':')
#                     # dis returns list of splitted text, [1] = name of dis, [2] = OMIM ID, extra chars need cleaning
#                     # print(dis[1][1:-5])
#                     # print(dis[2][:-1])
#                     id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})
#
#             # extract reactome ids and names
#             if 'Reactome;' in line:
#                 ract = line.split(';')
#                 # ract returns list with reactome id and name, needs cleaning
#                 id['Reactome'].update({ract[2][1:-2]: ract[1][1:]})
#                 # print(ract[1][1:])
#                 # print(ract[2][1:-2])
#
#             # look for functions
#             if ' F:' in line:
#                 if j < 5:
#                     # take only first 5 entries for now
#                     # print(j)
#                     fn = line.split(';')
#                     # fn returns list with GO ids and names
#                     id['Function'].update({fn[2][3:]: fn[1][1:]})
#                     # print(fn[1][1:])
#                     # print(fn[2][3:])
#                     j += 1
#
#             # look for biological processes
#             if ' P:' in line:
#                 if i < 5:
#                     # take only first 5 entries for now
#                     # print(i)
#                     bp = line.split(';')
#                     # bp returns list with GO ids and names
#                     id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
#                     #print(bp[1][1:])
#                     #print(bp[2][3:])
#                     i += 1
#
#             if 'GN   Name' in line:
#                 #print(line)
#
#                 if k == 0:
#                     gene = line.split('=')
#                     print(gene)
#                     gene = gene[1].split(' ')
#                     #print(gene[0])
#                     if ';' in gene[0]:
#                         gene=gene[0].split(';')
#                         #id['Gene'].update({gene[0]})
#                         gene = {'Gene': gene[0]}
#                         print(gene)
#                         #id.update(gene)
#                     else:
#                         gene = {'Gene':gene[0]}
#                         print(gene)
#                         #id.update(gene)
#                     id.update(gene)
#                     #print(id['Gene'])
#                     k +=1
#
#             #else:
#                 #id.update({'Gene': ''})
#
#             #print(id['Gene'])
#
#         Uniprot_Dict.append(id)
#
#     Uniprot_Dict = dict(zip(uniprot_id, Uniprot_Dict))
#
#     return(Uniprot_Dict)
    
def uniprot_rel(named_uprotList,itmpGraph):
    
    for item in named_uprotList:
        #print(named_uprotList[item]['Function'].keys())
        fun=list(named_uprotList[item]['Function'].keys())
        bp = list(named_uprotList[item]['BioProcess'].keys())
        for f in fun:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'],dict) :
                itmpGraph.add_association(Protein(namespace='HP',name=named_uprotList[item]['Gene']),BiologicalProcess(namespace='GOMF',name=f),
                                         citation='UniProt database',evidence='UniProt query')
            else:
                itmpGraph.add_association(Protein(namespace='HP',name=item),BiologicalProcess(namespace='GOMF',name=f),
                                         citation='UniProt database',evidence='UniProt query')
                
        
        for b in bp:
            if str(named_uprotList[item]['Gene']) != 'nan' and not isinstance(named_uprotList[item]['Gene'],dict):
                itmpGraph.add_association(Protein(namespace='HP',name=named_uprotList[item]['Gene']),BiologicalProcess(namespace='GOBP',name=b),
                                         citation='UniProt database',evidence='UniProt query')
            else:
                itmpGraph.add_association(Protein(namespace='HP',name=item),BiologicalProcess(namespace='GOBP',name=b),
                                         citation='UniProt database',evidence='UniProt query')
        
    return(itmpGraph)


def chem2gene2path_rel(named_chem2geneList, itmpGraph):
    for item in named_chem2geneList:
        # print(item)
        itemLen = len(named_chem2geneList[item]) - 1
        # print(itemLen)
        for j in range(itemLen - 1):
            # print(j)
            # print('yo',named_chem2geneList)
            # if 'accession' not in named_chem2geneList[j]:
            # continue

            itmpGraph.add_association(
                Protein(namespace='HP', name=named_chem2geneList[item][itemLen]['component_synonym']),
                BiologicalProcess(namespace='Reactome', name=named_chem2geneList[item][j]['xref_name']),
                citation='ChEMBL database', evidence='ChEMBL query',
                Reactome=named_chem2geneList[item][j]['xref_id'])

    return (itmpGraph)


def ExtractFromUniProt(uniprot_id):
    # from bioservices import UniProt
    Uniprot_Dict = []
    # Make a link to the UniProt webservice
    # service = UniProt()

    # uniprot_id = ['Q9Y5L3','P13716','Q89226299']
    mapped_uprot = []

    for id in uniprot_id:
        # print(id)

        # create URL for each uniprot id
        url = 'https://www.uniprot.org/uniprot/' + id + '.txt'
        print(url)

        #     #Retrieve data for id in text format if found in uniprot

        ret_uprot = requests.get(url)
        ret_uprot = ret_uprot.text.split('\n')
        # print(ret_uprot)

        # if len(ret_uprot) < 2:
        if ret_uprot == ['']:
            # uniprot_id.remove(id)
            # print(id)
            continue

        id_copy = id
        mapped_uprot.append(id_copy)
        i = 0
        j = 0
        k = 0
        id = {}
        id['Disease'] = {}
        id['Reactome'] = {}
        id['Function'] = {}
        id['BioProcess'] = {}
        id['Gene'] = {}
        # id['Gene'] = {}
        # print(id)

        # parse each line looking for info about disease, pathway, funcn, bp and so on
        for line in ret_uprot:
            # print(line)

            # line = line.decode('utf-8')

            # parse lines with disease and extract disease names and omim ids
            if '-!- DISEASE:' in line:
                if ('[MIM:' in line):
                    dis = line.split(':')
                    # dis returns list of splitted text, [1] = name of dis, [2] = OMIM ID, extra chars need cleaning
                    # print(dis[1][1:-5])
                    # print(dis[2][:-1])
                    id['Disease'].update({dis[1][1:-5]: dis[2][:-1]})

            # extract reactome ids and names
            if 'Reactome;' in line:
                ract = line.split(';')
                # ract returns list with reactome id and name, needs cleaning
                id['Reactome'].update({ract[2][1:-2]: ract[1][1:]})
                # print(ract[1][1:])
                # print(ract[2][1:-2])

            # look for functions
            if ' F:' in line:
                if j < 5:
                    # take only first 5 entries for now
                    # print(j)
                    fn = line.split(';')
                    # fn returns list with GO ids and names
                    id['Function'].update({fn[2][3:]: fn[1][1:]})
                    # print(fn[1][1:])
                    # print(fn[2][3:])
                    j += 1

            # look for biological processes
            if ' P:' in line:
                if i < 5:
                    # take only first 5 entries for now
                    # print(i)
                    bp = line.split(';')
                    # bp returns list with GO ids and names
                    id['BioProcess'].update({bp[2][3:]: bp[1][1:]})
                    # print(bp[1][1:])
                    # print(bp[2][3:])
                    i += 1

            if 'GN   Name' in line:
                # print(line)

                if k == 0:
                    gene = line.split('=')
                    # print(gene)
                    gene = gene[1].split(' ')
                    # print(gene[0])
                    if ';' in gene[0]:
                        gene = gene[0].split(';')
                        # id['Gene'].update({gene[0]})
                        gene = {'Gene': gene[0]}
                        print(gene)
                        # id.update(gene)
                    else:
                        gene = {'Gene': gene[0]}
                        # print(gene)
                        # id.update(gene)
                    id.update(gene)
                    # print(id['Gene'])
                    k += 1

                # else:
                # id.update({'Gene': ''})

                # print(id['Gene'])

        Uniprot_Dict.append(id)

    Uniprot_Dict = dict(zip(mapped_uprot, Uniprot_Dict))
    # len(uniprot_id)
    # print(Uniprot_Dict)

    return (Uniprot_Dict)