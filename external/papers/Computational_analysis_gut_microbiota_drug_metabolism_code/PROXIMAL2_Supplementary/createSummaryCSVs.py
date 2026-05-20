import pandas as pd
import time
from rdkit import Chem
from tqdm import tqdm
import json
from glob import glob

# Change which output directory to parse here
OUTPUT_DIRECTORY = 'output/set1_subset2_products'

def smiles2inchikey(smiles):
    mol = Chem.MolFromSmiles(smiles)
    key = Chem.inchi.MolToInchiKey(mol)
    return key

start_time = time.time()

# Create summary csvs    
print('Creating summary csvs...')

for drug_folder in tqdm(glob(f'{OUTPUT_DIRECTORY}/Met*')):
    drugid = drug_folder.split('_')[-1]
    drug_summary_df = None
    drug_names = []
    drug_smiles = []
    ecs = []
    reaction_centers = []
    product_smiles = []
    product_keys = []
    for reaction_json in glob(f"{drug_folder}/*/*"):
        reaction_data = json.load(open(reaction_json))
        for ec in reaction_data['TemplateReaction'][0]['ec'][0] + ['NoEC'] * ({0:1}.get(len(reaction_data['TemplateReaction'][0]['ec'][0]),0)):
            for rc in reaction_data['QueryInformation'][0]['reactionCenter(s)'] + ['NoRC'] * ({0:1}.get(len(reaction_data['QueryInformation'][0]['reactionCenter(s)']),0)):
                drug_names.append(reaction_data['QueryInformation'][0]['name'])
                drug_smiles.append(reaction_data['QueryInformation'][0]['smiles'])
                ecs.append(ec)
                reaction_centers.append(rc)
                product_smiles.append(reaction_data['GeneratedProduct'][0]['smiles'])
                product_keys.append(smiles2inchikey(reaction_data['GeneratedProduct'][0]['smiles']))
    drug_summary_df = pd.DataFrame({'DrugName': drug_names, 'DrugSmiles': drug_smiles, 'EC': ecs, 'ReactionCenter': reaction_centers, 'ProductSmiles': product_smiles, 'ProductInChiKey': product_keys})
    drug_summary_df.drop_duplicates().reset_index(drop=True).to_csv(f"{drug_folder}/{drugid}_output_summary.csv")

print('Done!')

print("--- %s seconds = %s hours ---" % ((time.time() - start_time),(time.time() - start_time)/60/60))
