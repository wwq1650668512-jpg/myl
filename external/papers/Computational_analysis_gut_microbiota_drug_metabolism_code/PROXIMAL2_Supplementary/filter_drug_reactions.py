import pandas as pd
import json

DATADIR = '../../drugbank_data' # directory with drugbank files

# import reactions parsed from xml download obtained from drugbank
rxns = pd.read_csv(f"{DATADIR}/drugbank_raw_reactions.csv") 

# And start by filtering out any reactions which explicitly mentione cytochrome or liver as well as reactions that do not have enzyme information
wEnzyme = rxns.dropna(subset=['Enzyme'])
nonLiver = wEnzyme.loc[wEnzyme.Enzyme.apply(lambda e: e != '' and  "cytochrome" not in e.lower() and "liver" not in e.lower())]

# import map of drugbankIDs to inchiKeys necessary for compound (both substrate and product) tractability
drugbank_productID2InchiKey = json.load(open(f'{DATADIR}/drugbank_productID2InchiKey.json'))

# we can only keep molecules for which we have information
prodInfo = nonLiver.loc[nonLiver.apply(lambda row: row['SubstrateID'] in drugbank_productID2InchiKey and row['ProductID'] in drugbank_productID2InchiKey, axis=1)]

uniLinks = pd.read_csv('../linkFiles/uniprot_links.csv') # map of names of enzymes to uniProt IDs (provided in github)
keggIds = json.load(open('../linkFiles/uniprotID_to_KEGGECs_v2.json')) # map of uniProt to EC numbers (provided in github)

keggDict = {}
for d in keggIds:
    keggDict[d['UniprotID']] = d['EC']


# Last level of filtering is to filter for any reactions where we cannot obtain the EC number associated with the enzyme
# as we need that for tractability of enzymes

subList = []
prodList = []
ecList = []

for idx, row in prodInfo.iterrows():
    possibleUniProtIds = uniLinks.loc[uniLinks['UniProt Name'] == row['Enzyme']]['UniProt ID']
    for upID in possibleUniProtIds:
        possibleECs = keggDict.get(upID, [])
        for ec in possibleECs:
            subList.append(row['SubstrateID'])
            prodList.append(row['ProductID'])
            ecList.append(ec)

final_rxns = pd.DataFrame({'SubstrateID': subList, 'ProductID': prodList, 'Enzyme': ecList}).drop_duplicates()

print(f"After filtering the E-Biological set contains {len(final_rxns)} reactions")



