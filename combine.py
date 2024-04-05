import os
import json

def combine_json_files(directory, output_file):
    combined_data = {}
    for filename in os.listdir(directory):
        if filename.endswith('.json'):
            file_path = os.path.join(directory, filename)
            with open(file_path, 'r') as file:
                data = json.load(file)
                combined_data.update(data)
    with open(output_file, 'w') as output:
        json.dump(combined_data, output,indent=4, sort_keys=True)

# Example usage:
directory = '/home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json/UniProtKB_separate'
output_file = '/home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json/UniProtKB_dict.json'

combine_json_files(directory, output_file)
