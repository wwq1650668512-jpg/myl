import json
import csv
from os.path import isfile

import torch
from torch_geometric.data import Data
import requests

from gnn_som import createGnnSom, loadGnnSomState
from gnn_som.MolFromKcf import MolFromKcfFile

def config_features(mol, enzyme, config):
    numFeatures = sum(len(feature) for feature in config['features'].values())
    x = torch.zeros((mol.GetNumAtoms(), numFeatures), dtype=torch.float32)
    for atom in mol.GetAtoms():
        x[atom.GetIdx(), config['features']['enzyme'].index(enzyme)] = 1
        offset = len(config['features']['enzyme'])
        x[atom.GetIdx(), offset + config['features']['element'].index(atom.GetSymbol())] = 1
        offset += len(config['features']['element'])
        x[atom.GetIdx(), offset + config['features']['kcfType'].index(atom.GetProp('kcfType'))] = 1
    
    edgeIndex = torch.zeros((2, mol.GetNumBonds() * 2), dtype=torch.int64)
    for bond in mol.GetBonds():
        i = bond.GetIdx()
        edgeIndex[0][i * 2] = bond.GetBeginAtomIdx()
        edgeIndex[1][i * 2] = bond.GetEndAtomIdx()
        edgeIndex[0][i * 2 + 1] = bond.GetEndAtomIdx()
        edgeIndex[1][i * 2 + 1] = bond.GetBeginAtomIdx()
    
    return x, edgeIndex

def get_likelihood_scores(conversionFile, inputFile):
    drug_smiles_list = []
    ec_full_list = []
    ec_first_list = []
    ec_second_list = []
    reaction_centers_list = []
    drug = []
    predicted_product = []
    product_inchi = []
    products_library = [] #library to store SOM value of each product
    key_value_pairs = {}
    
    # Store KEGG ID of each drug
    with open(conversionFile, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            key = row['name']
            value = row['ID']
            key_value_pairs[key] = value
            
    with open(inputFile, 'r') as f:
        data = json.load(f)
        for key, value in data.items():
            underscore_index = key.rfind("_")
            drug_name=(key[:underscore_index])
            
            #Ipsalazide is not in KEGG library
            if not (drug_name == "Ipsalazide" or "Neoprontosil"):
                if drug_name in key_value_pairs:
                    kegg_id = key_value_pairs[drug_name]
                    drug.append(kegg_id)
    
                predicted_product.append(key[underscore_index + 1:])
                drug_smiles_list.append(value["Drug_Smiles"])
                ec_full_list.append(value["EC_Full"])
                ec_first_list.append(value["EC_First"])
                ec_second_list.append(value["EC_Second"])
                product_inchi.append(value["Product_InChi_Key"])
                reaction_centers_list.append(value["Reaction_Centers"])
            
    # Configurate 10 GNN-SOM models to apply
    with open('data/config.json', 'r') as f:
        config = json.load(f)
    config['features']['enzyme'] = [tuple(ec) for ec in config['features']['enzyme']] 
    
    models = []
    for i, params in enumerate(config['models']):
        model = createGnnSom(*config['models'][i])
        loadGnnSomState(model, torch.load('data/model%d.pt' % i, map_location=torch.device('cpu')))
        models.append(model)
    

    for a in range(len(drug)):
        molecule = drug[a]
        enzyme = (ec_first_list[a], ec_second_list[a])
    
        if not isfile(molecule + '.kcf'):
            with open(molecule + '.kcf', 'wb') as f:
                f.write(requests.get('https://www.genome.jp/entry/-f+k+' + molecule).content)
        mol = MolFromKcfFile(molecule + '.kcf')  
        
        x, edgeIndex = config_features(mol, enzyme, config)
        data = Data(x=x, edgeIndex=edgeIndex)

        y = None
        for model in models:
            newY = torch.sigmoid(model(data.x, data.edgeIndex))
            y = newY if y is None else torch.add(y, newY)
        y = torch.div(y, len(models))
      
        SOM_value = 0
        for atom in mol.GetAtoms():
            label = y[atom.GetIdx()].item()
            #Take max SOM value if there are multiple reaction centers 
            for z in range(len(reaction_centers_list[a])):
                if (atom.GetIdx() == reaction_centers_list[a][z]):
                    if (float(label) > SOM_value):
                        SOM_value = float(label)

        dict = {
            "Drug": drug[a],
            "Product": predicted_product[a],
            "EC": ec_full_list[a],
            "Value": SOM_value,
            "InChi": product_inchi[a]}
        products_library.append(dict)
        
    return products_library
