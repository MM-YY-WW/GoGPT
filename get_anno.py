import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
import pandas as pd
import csv
import os
from tqdm import tqdm
def get_sorted_unique_goall(gene_id):
    robjects.r('library(org.Hs.eg.db)')
    robjects.r('library(base)')
    r_code = f'sort(unique(select(org.Hs.eg.db, keys = c("{gene_id}"), columns = c("GOALL"), keytype = "ENTREZID")$GOALL))'
    sorted_unique_goall = robjects.r(r_code)
    sorted_unique_goall = list(sorted_unique_goall)
    return sorted_unique_goall

def parse_obo_file(file_path, go_id_list):
    go_terms = {}
    current_term = None

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()

            if line.startswith('[Term]'):
                if current_term is not None and current_term['id'] in go_id_list:
                    go_terms[current_term['id']] = current_term
                current_term = {'id': None, 'name': None, 'namespace': None, 'def': None, 'synonym': []}
            elif line.startswith('id:'):
                current_term['id'] = line.split(': ')[1]
            elif line.startswith('name:'):
                current_term['name'] = line.split(': ')[1]
            elif line.startswith('namespace:'):
                current_term['namespace'] = line.split(': ')[1]
            elif line.startswith('def:'):
                current_term['def'] = line.split('"')[1]
            elif line.startswith('synonym:'):
                current_term['synonym'].append(line.split('"')[1])

    if current_term is not None and current_term['id'] in go_id_list:
        go_terms[current_term['id']] = current_term

    return go_terms

def create_paragraph(go_terms):
    paragraph = ""
    for go_id, term in go_terms.items():
        paragraph += term["name"] + ". "  # Add the Name of the term to the paragraph
    return paragraph.strip()

def save_unique_gene_ids(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
    unique_gene_ids = sorted(df['NCBI GeneID'].unique())
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['NCBI GeneID'])
        for gene_id in unique_gene_ids:
            writer.writerow([gene_id])

def save_gene_id_paragraphs(input_file, output_file, obo_file_path):
    gene_ids = pd.read_csv(input_file)['NCBI GeneID']
    for gene_id in tqdm(gene_ids[50000:60000]):
        paragraph = create_paragraph(parse_obo_file(obo_file_path, get_sorted_unique_goall(str(gene_id))))
        if paragraph != "":
            df = pd.DataFrame({'NCBI GeneID': [gene_id], 'Paragraph': [paragraph]})
            df.to_csv(output_file, index=False, mode='a', header=not os.path.exists(output_file)) 
            
obo_file_path = '/home/yuwei/codes/Gene_ontology_PT/go-basic.obo' 
tsv_file_path = '/home/yuwei/codes/Gene_ontology_PT/ncbi_dataset.tsv'
csv_unique_gene_ids = '/home/yuwei/codes/Gene_ontology_PT/unique_gene_ids.csv'
csv_gene_id_paragraphs = '/home/yuwei/codes/Gene_ontology_PT/entrezid_gotext50000_60000.csv'

save_unique_gene_ids(tsv_file_path, csv_unique_gene_ids)

save_gene_id_paragraphs(csv_unique_gene_ids, csv_gene_id_paragraphs, obo_file_path)

