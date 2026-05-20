import json
from itertools import groupby
import csv

# Initialize an empty dictionary to store key-value pairs
key_value_pairs = {}

# Open the CSV file and read its content
def drugID_name_dictionary(nameConversionFile):
    with open(nameConversionFile, 'r') as csv_file:
        # Create a CSV reader object
        csv_reader = csv.DictReader(csv_file)
    
        # Iterate through each row in the CSV file
        for row in csv_reader:
            key = row['name']
            value = row['ID']
            key_value_pairs[key] = value
    return key_value_pairs

# Create a library storing each drugs list of inchis and their ids to later comb 
# through for identical products/inchis
def create_product_library(inputFile):
    sublib = list()
    with open(inputFile, 'r') as f:
        all_drugs = list()
        data = json.load(f)
        previous = "starter"
        drug_dict = {
            'Name' : "Ignore",
            'Products_ID' : list(),
            'Products' : list()
            }
        
        for key,value in data.items():
            underscore_index = key.rfind("_")
            drug = (key[:underscore_index])
            index = key[underscore_index+1:]
            product = value['Product_InChi_Key'][:14]
            if (not drug == previous):
                sublib.append(drug_dict)
                drug_dict = {
                    'Name' : drug,
                    'Products_ID' : list(),
                    'Products' : list()
                    }
            drug_dict['Products'].append(product)
            drug_dict['Products_ID'].append(index)
            previous = drug
        del sublib[0] #Remove placeholder value
    return sublib
        
# Function to identify position of identical product InChis
def find_identical_positions(lst):
    positions = {}
    
    for i, item in enumerate(lst):
        if item not in positions:
            positions[item] = [i]
        else:
            positions[item].append(i)
    
    # Filter positions for identical strings
    identical_positions = {key: value for key, value in positions.items() if len(value) > 1}
    return identical_positions


    
def create_unique_product_library(sublib):
    # Create final library with list of ids of identical positions for each drug
    new_library = list()

    for drug in sublib:
        drug_dict = {
            'Name' : drug['Name'],
            'Identical Positions': list()
            }
        identical = find_identical_positions((drug['Products']))
        correct_id = drug['Products_ID']
            
        for product in identical:
            position_list = identical[product]
            for i in range (len(position_list)):
                position = position_list[i]
                position_list[i] = correct_id[position]
            drug_dict['Identical Positions'].append(position_list)
        new_library.append(drug_dict)
    return new_library
    

def remove_duplicate_products(inputFile, nameConversionFile, group):
    key_value_pairs = drugID_name_dictionary(nameConversionFile)
    sublib = create_product_library(inputFile)
    new_library = create_unique_product_library(sublib)
    
    id_list = list()
    drug = group[0]['Drug']
    for entry in new_library:
        # kegg_id = entry["Name"]
        kegg_id = key_value_pairs[entry["Name"]]
        if kegg_id == drug:
            id_list = entry["Identical Positions"]
    for duplicate_list in id_list:
        max_score = 0
        product_positions = list()
        for position in duplicate_list:
            #Get the score of it and update max or remove it
            for i in range (len(group)):
                drug = group[i]
                if drug["Product"] == position:
                    score = drug['Value']
                    if (score > max_score):
                        max_score = score
                    product_positions.append(i)
        
        # Note positions of all products that were not the highest score
        delete_this = list()
        for j in range (len(product_positions)):
            pos = product_positions[j]
            if (not group[pos]['Value'] == max_score):
                # print(max_score, "and", group[j]['Value'])
                delete_this.append(pos)
        
        delete_this = sorted(delete_this, reverse=True)
        for k in range (len(delete_this)):
            position = delete_this[k]
            # print(group[position])
            del group[position]
    return group

def rank_GNN_Results(valuesFile, outputFile, inputFile, nameConversionFile):
    with open(valuesFile, 'r') as f:
        data = json.load(f)
       # Sort the data by Drug and then by Value within each Drug group
        sorted_data = sorted(data, key=lambda x: (x['Drug'], x['Value']), reverse=True)
        
        # Group the sorted data by Drug
        grouped_data = {key: list(group) for key, group in groupby(sorted_data, key=lambda x: x['Drug'])}
        
        # Rank the grouped data by Value
        # pos = 1
        with open(outputFile, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            
            header = ['Drug', 'Product ID', 'EC', 'Raw Score', 'Rank']  
            csv_writer.writerow(header)
            
            for drug, group in grouped_data.items():
                group = remove_duplicate_products(inputFile, nameConversionFile, group)
                total = len(group)
                ranked_group = sorted(group, key=lambda x: x['Value'], reverse=True)
                print(f"Group: {drug}")
        
                for rank, item in enumerate(ranked_group, start=1):
                    if (rank >= 1):
                           column1 = item['Drug']
                           column2 = item['Product']
                           column3 = item['EC']
                           column4 = item['Value']
                           column5 = str(rank) + "/"  + str(total)
                           csv_writer.writerow([column1, column2, column3, column4, column5])
                           print(f"  Rank {rank}/{total}: Drug: {item['Drug']} - Product ID: {item['Product']} - Value: {item['Value']}")
            
inputFile = 'data/updated_set2_gnn_input.json'          
valuesFile = 'data/GNN_values_set2_upload2.json'
outputFile = 'data/set2Ranked_upload.csv'
nameConversionFile = 'data/Downloads/conv2.csv'  

rank_GNN_Results(valuesFile, outputFile, inputFile, nameConversionFile)