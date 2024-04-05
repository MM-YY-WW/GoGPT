import pandas as pd
import requests
import os
from tqdm import tqdm
from lxml import html
import urllib.request
import json
from bs4 import BeautifulSoup as BS
import urllib

id_dict= {}
def load_json(DB_list):
    folder_path= '/home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json'
    full_dict = {}
    for i in DB_list:
        json_path = os.path.join(folder_path, i+'_dict.json')
        if os.path.exists(json_path):
            with open(json_path, 'r') as file:
                full_dict[i] = json.load(file)
                print(f"Load {i} Dictionary")
        else:
            full_dict[i] = {}
            print(f"Create {i} Dictionary")

    return full_dict

def save_json(full_dict, DB_list, filename):
    #folder_path= '//home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json/UniProtKB_separate'
    folder_path= '//home/yuwei/codes/Gene_ontology_PT/ID_Seq_Json'
    for i in DB_list:
        #json_path = os.path.join(folder_path, i+'_'+filename[:4]+'_dict.json')
        json_path = os.path.join(folder_path, i+'_dict.json')
        with open(json_path, 'w') as file:
            json.dump(full_dict[i], file, indent=4, sort_keys=True)
            print(f"Save {i} Dictionary")

# dict{'DB_Object_ID':[DNA_seq, RNA_seq, Pr_seq]}

def iter_folder(foldername):
    for filename in os.listdir(foldername):
        file_path = os.path.join(foldername, filename)

        #if os.path.isfile(file_path) and filename =='genedb_lmajor.gaf' and  filename !='ecocyc.gaf' and filename != 'cgd.gaf':
        if os.path.isfile(file_path) and filename == 'xenbase.gaf':
        #if os.path.isfile(file_path):
            iter_gaf_file(file_path, filename)
        

def iter_gaf_file(filepath, filename):
    gaf_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'Qualifier', 'GO_ID', 'DB:Reference', 'Evidence_Code', 'With_or_From', 'Aspect', 'DB_Object_Name', 'DB_Object_Synonym', 'DB_Object_Type', 'Taxon', 'Date', 'Assigned_By', 'Annotation_Extension', 'Gene_Product_Form_ID']
    gaf_data = pd.read_csv(filepath, sep='\t', comment='!', names=gaf_columns)
    outfilename = os.path.join('/home/yuwei/codes/Gene_ontology_PT/Seq_Gotext',filename[:-4]+'.csv')
    output_columns = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'DNA_seq', 'RNA_seq', 'Pr_seq', 'GO_ID']
    outputfile=pd.DataFrame(columns=output_columns)
    #print(filename,len(gaf_data))
    #print(filename,list(set(gaf_data['DB'])))
    full_dict=load_json(list(set(gaf_data['DB'])))

    #resume_from_index = 207394
    resume_from_index = 0
    for i in tqdm(range(resume_from_index, len(gaf_data))):
    #for i in tqdm(range(10)):
        outputfile.loc[len(outputfile.index)], full_dict=get_sequence(gaf_data.iloc[i], full_dict)
        if len(outputfile) == 100 or i == len(gaf_data)-1 :
                outputfile.to_csv(outfilename, index=False, mode='a',header=not os.path.exists(outfilename))
                outputfile=pd.DataFrame(columns=output_columns)
                save_json(full_dict, list(set(gaf_data['DB'])), filename)

def get_sequence(gaf_obj, full_dict):
    db = gaf_obj['DB']
    db_object_id = gaf_obj['DB_Object_ID']
    db_object_symbol = gaf_obj['DB_Object_Symbol']
    Go_id = gaf_obj['GO_ID']
    #print(gaf_obj)
    if db == 'ComplexPortal':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_ComplexPortal(db_object_id, db_object_symbol, full_dict[db])
    elif db == 'CGD':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_CGD(db_object_id, db_object_symbol,full_dict[db])
    elif db == 'PomBase':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_PomBase(db_object_id, db_object_symbol,full_dict[db])
    elif db == 'dictyBase':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_dictyBase(db_object_id, db_object_symbol,full_dict[db])
    elif db == 'UniProtKB':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_UniProtKB(db_object_id, db_object_symbol,full_dict[db])
    elif db == 'TriTrypDB':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_TritrypDB(db_object_id, db_object_symbol,full_dict[db])
    elif db == 'Xenbase':
        DNA_seq, RNA_seq, Pr_seq, dict = retrieve_Xenbase(db_object_id, db_object_symbol,full_dict[db])
        
    else:
        raise NotImplementedError(f"{db} is not Implemented !")

    full_dict[db] = dict
    return [db, db_object_id, db_object_symbol, DNA_seq, RNA_seq, Pr_seq, Go_id], full_dict

def retrieve_ComplexPortal(DB_Object_ID, DB_Object_Symbol,dict):
    '''
    For dataset:
    goa_cow_complex.gaf

    '''
    url = f'https://www.ebi.ac.uk/intact/complex-ws/export/{DB_Object_ID}'
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if DB_Object_ID not in dict.keys():
        dict[DB_Object_ID] = ['','','']
        try:
            response = requests.get(url)
            for i in range(len(response.json()['data'])):
                if 'sequence' in response.json()['data'][i].keys():
                    Pr_seq = response.json()['data'][i]['sequence']
            print(DB_Object_ID)
            print(Pr_seq)
            dict[DB_Object_ID][2] = Pr_seq
        except Exception as e:
            print(e)
    else:
        Pr_seq = dict[DB_Object_ID][2]
    return DNA_seq, RNA_seq, Pr_seq, dict

def retrieve_CGD(DB_Object_ID, DB_Object_Symbol, dict):
    '''
    For dataset:
    cgd.gaf (363843)
    '''
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if 'CAWG' in DB_Object_ID:
        print(DB_Object_ID)
        if DB_Object_ID not in id_dict.keys():
            CAWG_url =f'http://www.candidagenome.org/cgi-bin/locus.pl?locus={DB_Object_ID}&organism=C_albicans_SC5314'
            try:
                CAWG_response = requests.get(CAWG_url)
                for i in CAWG_response.text.split('http://www.candidagenome.org/cgi-bin/reference/litGuide.pl?dbid='):
                    if i[:3] == 'CAL':
                        id_dict[DB_Object_ID] =  i[:13]
                        DB_Object_ID = i[:13]
                        print(DB_Object_ID)
                        break
                if 'CAWG' in DB_Object_ID:
                    return DNA_seq, RNA_seq, Pr_seq, dict
            except Exception as e:
                print(e)
        else:
            DB_Object_ID = id_dict[DB_Object_ID]
    elif 'CORT' in DB_Object_ID or 'DEHA' in DB_Object_ID or'LELG' in DB_Object_ID or 'CTRG' in DB_Object_ID or 'CLUG' in DB_Object_ID or 'PGUG' in DB_Object_ID:
        print(DB_Object_ID)
        if DB_Object_ID not in id_dict.keys():
            SC5314_url = f'http://www.candidagenome.org/cgi-bin/search/textSearch?query={DB_Object_ID}&type=homolog&organism=C_albicans_SC5314'
            try:
                SC5314_response = requests.get(SC5314_url)
                soup=BS(SC5314_response.text, 'html.parser')
                ID_url = soup.find_all('a', href=lambda href:(href and href.startswith('http://www.candidagenome.org/cgi-bin/locus.pl?locus=')))[0]['href'] 
                try:
                    ID_response = requests.get(ID_url)
                    soup = BS(ID_response.content, 'html.parser')
                    for i in ID_response.text.split('http://www.candidagenome.org/cgi-bin/reference/litGuide.pl?dbid='):
                        if i[:3] == 'CAL':
                            id_dict[DB_Object_ID] =  i[:13]
                            DB_Object_ID = i[:13]
                            print(DB_Object_ID)
                            break
                    if 'CORT' in DB_Object_ID or 'DEHA' in DB_Object_ID or'LELG' in DB_Object_ID or 'CTRG' in DB_Object_ID or 'CLUG' in DB_Object_ID or 'PGUG' in DB_Object_ID:
                        return DNA_seq, RNA_seq, Pr_seq, dict
                except Exception as e:
                    print(e)
            except Exception as e:
                print(e)
        else:
            DB_Object_ID = id_dict[DB_Object_ID]
    elif 'OVF' in DB_Object_ID or 'CJI97' in DB_Object_ID:
        print(DB_Object_ID)
        return DNA_seq, RNA_seq, Pr_seq, dict
    
    if DB_Object_ID not in dict.keys():
        dict[DB_Object_ID] = ['','','']
        print(DB_Object_ID)

        DNA_url = f'http://www.candidagenome.org/cgi-bin/getSeq?map=a3map&seq={DB_Object_ID}'
        Pr_url=f'http://www.candidagenome.org/cgi-bin/getSeq?map=p3map&seq={DB_Object_ID}'
        try:
            DNA_response = requests.get(DNA_url)
            DNA_tree = html.fromstring(DNA_response.content)
            DNA_seq = ''.join(DNA_tree.xpath('/html/body/pre')[0].text.split('\n')[1:])
            dict[DB_Object_ID][0] = DNA_seq
            print(DNA_seq)
        except Exception as e:
            print(e)
        try:
            Pr_response = requests.get(Pr_url)
            Pr_tree = html.fromstring(Pr_response.content)
            Pr_seq = ''.join(Pr_tree.xpath('/html/body/pre')[0].text.split('\n')[1:])
            dict[DB_Object_ID][2] = Pr_seq
            print(Pr_seq)
        except Exception as e:
            print(e)
    else:
        DNA_seq = dict[DB_Object_ID][0]
        Pr_seq = dict[DB_Object_ID][2]
    return DNA_seq, RNA_seq, Pr_seq, dict

def retrieve_PomBase(DB_Object_ID, DB_Object_Symbol,dict):
    '''
    
    '''
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if DB_Object_ID not in dict.keys():
        dict[DB_Object_ID] = ['','','']
        Pr_seq = get_sequence_by_id(fasta_file='/home/yuwei/codes/Gene_ontology_PT/Downloaded_sequence_files/PomBase/peptide.fa', search_id=DB_Object_ID)
        DNA_seq = get_sequence_by_id(fasta_file='/home/yuwei/codes/Gene_ontology_PT/Downloaded_sequence_files/PomBase/cds+introns+utrs.fa', search_id=DB_Object_ID)
        dict[DB_Object_ID][0] = DNA_seq
        dict[DB_Object_ID][2] = Pr_seq
        print(DB_Object_ID)
        print(DNA_seq)
        print(Pr_seq)
    else:
        DNA_seq = dict[DB_Object_ID][0]
        Pr_seq = dict[DB_Object_ID][2]

    return DNA_seq, RNA_seq, Pr_seq, dict   

def read_fasta(filename):
    sequences = {}
    current_id = None
    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ''
            else:
                sequences[current_id] += line
    return sequences

def get_sequence_by_id(fasta_file, search_id):
    seq = ''
    sequences = read_fasta(fasta_file)
    for i in range(len(sequences.keys())):
        if search_id in list(sequences.keys())[i]:
            seq = sequences.get(list(sequences.keys())[i], '')
    return seq

def retrieve_dictyBase(DB_Object_ID, DB_Object_Symbol,dict):
    '''
    For dataset:
    dictbase.gdf(81452)
    '''
    url = f'http://dictybase.org/gene/{DB_Object_ID}'
    Pr_seq, DNA_seq, RNA_seq, id = '','','',''

    if DB_Object_ID not in dict.keys():
        try:
            response = requests.get(url)
            tree = html.fromstring(response.content)
            for i in tree.xpath("//script[contains(text(), 'var config')]")[0].text.split('/protein/'):
                if i[0:3]=='DDB':
                    id = i.split('.json')[0]
                    break
            print(DB_Object_ID)
            print(id)
            dict[DB_Object_ID] = ['','','']
        except Exception as e:
            print(e)
        Pr_url = f'http://dictybase.org/db/cgi-bin/dictyBase/yui/get_fasta.pl?decor=1&primary_id={id}&sequence=Protein'
        DNA_url = f'http://dictybase.org/db/cgi-bin/dictyBase/yui/get_fasta.pl?decor=1&primary_id={id}&sequence=DNA%20coding%20sequence'
        try:
            Pr_response = requests.get(Pr_url)
            Pr_seq = ''.join(Pr_response.text.split('\n')[1:-1])
            print(Pr_seq)
            dict[DB_Object_ID][2] = Pr_seq
        except Exception as e:
            print(e)
        try:
            DNA_response = requests.get(DNA_url)
            DNA_seq = ''.join(DNA_response.text.split('\n')[1:-1])
            print(DNA_seq)
            dict[DB_Object_ID][0] = DNA_seq
        except Exception as e:
            print(e)
    else:
        Pr_seq = dict[DB_Object_ID][2]
        DNA_seq = dict[DB_Object_ID][0]
    
    return DNA_seq, RNA_seq, Pr_seq, dict

def retrieve_UniProtKB(DB_Object_ID, DB_Object_Symbol, dict):
    '''
    For dataset:
    ecocyc.gaf
    '''
    url = f'https://rest.uniprot.org/uniprotkb/{DB_Object_ID}.fasta'
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if DB_Object_ID not in dict.keys():  
        dict[DB_Object_ID] = ['','',''] 
        print(DB_Object_ID)
        try:
            response = requests.get(url)
            Pr_seq = ''.join(response.text.split("\n")[1:])
            dict[DB_Object_ID][2] = Pr_seq
            print(Pr_seq)
        except Exception as e:
            print(e)
    else:
        Pr_seq = dict[DB_Object_ID][2]
    return DNA_seq, RNA_seq, Pr_seq, dict
    
def retrieve_TritrypDB(DB_Object_ID, DB_Object_Symbol, dict):
    '''
    https://tritrypdb.org/tritrypdb/app/downloads
    '''
    if 'LmjF' in DB_Object_ID:
        DNA_path='/home/yuwei/codes/Gene_ontology_PT/Downloaded_sequence_files/TritrypDB/TriTrypDB-24_LmajorFriedlin_AnnotatedCDSs.fasta'
        RNA_path='/home/yuwei/codes/Gene_ontology_PT/Downloaded_sequence_files/TritrypDB/TriTrypDB-24_LmajorFriedlin_AnnotatedTranscripts.fasta'
        Pr_path = '/home/yuwei/codes/Gene_ontology_PT/Downloaded_sequence_files/TritrypDB/TriTrypDB-24_LmajorFriedlin_AnnotatedProteins.fasta'
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if DB_Object_ID not in dict.keys(): 
        dict[DB_Object_ID] = ['','',''] 
        print(DB_Object_ID)
        DNA_seq = get_sequence_by_id(fasta_file=DNA_path, search_id=DB_Object_ID)
        print(DNA_seq)
        RNA_seq = get_sequence_by_id(fasta_file=RNA_path, search_id=DB_Object_ID)
        print(RNA_seq)
        Pr_seq = get_sequence_by_id(fasta_file=Pr_path, search_id=DB_Object_ID)
        print(Pr_seq)
        dict[DB_Object_ID][0],dict[DB_Object_ID][1],dict[DB_Object_ID][2]=DNA_seq, RNA_seq, Pr_seq
    else:
        DNA_seq, RNA_seq, Pr_seq = dict[DB_Object_ID][0],dict[DB_Object_ID][1],dict[DB_Object_ID][2]
    return DNA_seq, RNA_seq, Pr_seq, dict

def retrieve_Xenbase(DB_Object_ID, DB_Object_Symbol, dict):
    Pr_seq, DNA_seq, RNA_seq = '','',''
    if DB_Object_ID not in dict.keys(): 
        dict[DB_Object_ID] = ['','',''] 
        print(DB_Object_ID)
    
    
    else:
        DNA_seq, RNA_seq, Pr_seq = dict[DB_Object_ID][0],dict[DB_Object_ID][1],dict[DB_Object_ID][2]
    return DNA_seq, RNA_seq, Pr_seq, dict
iter_folder('/home/yuwei/codes/Gene_ontology_PT/Go_annotation')

#/html/body/div[4]