import csv
import json

columns = ['DB','DB_Object_ID','DB_Object_Symbol','DNA_seq','RNA_seq','Pr_seq','GO_ID']

csv_file_path = '/home/yuwei/codes/Gene_ontology_PT/Seq_Gotext/cgd.csv'
json_file_path = '/home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json/CGD_dict.json'

full_dict = {}

with open(csv_file_path, mode='r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        full_dict[row['DB_Object_ID']] = [row['DNA_seq'], row['RNA_seq'], row['Pr_seq']]

with open(json_file_path, 'w') as file:
    json.dump(full_dict, file, indent=4, sort_keys=True)
