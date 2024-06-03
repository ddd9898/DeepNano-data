import pandas as pd 
import numpy as np
import pymol
from pymol import *
from Bio import SeqIO
import os
from tqdm import tqdm

aa_dict = {'GLY':'G','ALA' :'A','VAL': 'V','LEU': 'L','ILE': 'I','PRO': 'P','PHE' :'F','TYR' :'Y','TRP' :'W',
'SER' :'S','THR': 'T','CYS': 'C','MET': 'M','ASN': 'N','GLN' :'Q','ASP': 'D', 'GLU' :'E','LYS' :'K','ARG' :'R', 'HIS' :'H',
'NAG':'X'}

# antigen_dict = dict()
# antigen_Emd_dict = dict()
# antigen_data = list()
# for fa in SeqIO.parse('./data/nanobody_20221115/unique_antigen.fasta',  "fasta"):
#     seqs = str(fa.seq).split('\t')[0]
#     ID = str(fa.description)
#     antigen_data.append([ID,seqs])
#     antigen_dict[ID] = seqs



dis = 5

metaData = pd.read_csv('./all_nano_structures/sabdab_nano_summary_all.tsv',sep='\t',keep_default_na=False).values.tolist()
# metaData = pd.read_csv('./nanobody_pan3/data/Sabdab/all_nano_structures/sabdab_nano_summary_all.tsv',sep='\t',keep_default_na=False).values.tolist()

seq_list = list()
unique_seq = list()
for item in tqdm(metaData):
    pdb = item[0]
    Hchain = item[1]
    Lchain = item[2]
    antigen_chain = item[4]
    antigen_type = item[5]
    antigen_name = item[7]
    affinity = item[25]
    affinity_method	= item[27]
    
    if (Lchain != 'NA') or (antigen_chain == Hchain):
        continue
    
    if antigen_chain == 'NA' or '|' in antigen_chain:
        continue
    
    cmd.delete('all')
    cmd.load('./all_nano_structures/{}.pdb'.format(pdb))
    cmd.remove('solvent')
    cmd.remove('hetatm')
    
    #######################Step 1:Get nanobody Heavy chain#######################
    cmd.save('./{}.fasta'.format(pdb),'chain {}'.format(Hchain))
    for fa in SeqIO.parse('./{}.fasta'.format(pdb),  "fasta"):
        seq_Hchain = str(fa.seq)
    os.remove('./{}.fasta'.format(pdb))
    if seq_Hchain.count('?')>3:
        cmd.delete('all')
        continue
    else:
        seq_Hchain = seq_Hchain.replace('?','X')
        
    #Delete items with len(seq_nanobody) > 150
    if len(seq_Hchain)>150:
        continue
    
    
    #Get nanobody residues within dis of antigen (nanobody_pockets)
    seq_nanobody_pockets = list()
    
    #Get all residues indexs on nanobody_chain
    resi_index_list_all = {'list':[]}
    cmd.iterate('chain {}'.format(Hchain),"list.append((resi,resn))",space = resi_index_list_all)
    resi_index_list_temp = [item for item in list(set(resi_index_list_all['list']))]
    resi_index_list_temp.sort(key = resi_index_list_all['list'].index)

    resi_index_list_all = list()
    seq_nanobody_new = list()
    count = 0
    for m,item in enumerate(resi_index_list_temp):
        if item[1] not in aa_dict.keys():
            if 'X' == seq_Hchain[count]:
                seq_nanobody_new.append('X')
                resi_index_list_all.append(item[0])
                count += 1
        else:
            if aa_dict[item[1]] == seq_Hchain[count]:
                seq_nanobody_new.append(aa_dict[item[1]])
                resi_index_list_all.append(item[0])
                count += 1
        if count >= len(seq_Hchain):
            break
    seq_nanobody_new = ''.join(seq_nanobody_new)
    
    if seq_nanobody_new != seq_Hchain:
        print("Something is wrong with nanobody sequence")
        cmd.delete('all')
        continue
    
    #Select pocket indexs
    cmd.select('temp','byres chain {} within {} of chain {}'.format(Hchain,dis,antigen_chain))
    cmd.save('./temp.fasta','temp')
    for fa in SeqIO.parse('./temp.fasta',  "fasta"):
        seq_nanobody_pockets.append(Hchain+':'+str(fa.seq))
    os.remove('./temp.fasta')
    
    resi_index_list = {'list':[]}
    cmd.iterate('temp',"list.append((resi,resn))",space = resi_index_list)
    resi_index_list_temp = [item for item in list(set(resi_index_list['list']))]
    resi_index_list_temp.sort(key = resi_index_list['list'].index)
    resi_index_list = ','.join([str(resi_index_list_all.index(item[0])) for item in resi_index_list_temp if item[0] in resi_index_list_all])
    # resi_index_list = [resi_index_list_all.index(item[0]) for item in resi_index_list_temp if item[0] in resi_index_list_all]
    # seq_antigen_pockets_new = ''.join([seq_antigen[item] for item in resi_index_list])
    binding_site_nanobody = resi_index_list
    
    
    #####################Step 2: Get antigen chain########################
    cmd.save('./temp.fasta','chain {}'.format(antigen_chain))
    for fa in SeqIO.parse('./temp.fasta',  "fasta"):
        seq_antigen = str(fa.seq)
    os.remove('./temp.fasta')
    LEN_SEQ_ANTIGEN = len(seq_antigen)
    if seq_antigen.count('?')>3:
        cmd.delete('all')
        continue
    else:
        seq_antigen = seq_antigen.replace('?','X')
    

    #Get antigen residues within dis of nanobody (antigen_pockets)
    seq_antigen_pockets = list()
    
    #Get all residues indexs on antigen_chain
    resi_index_list_all = {'list':[]}
    cmd.iterate('chain {}'.format(antigen_chain),"list.append((resi,resn))",space = resi_index_list_all)
    resi_index_list_temp = [item for item in list(set(resi_index_list_all['list']))]
    resi_index_list_temp.sort(key = resi_index_list_all['list'].index)

    resi_index_list_all = list()
    seq_antigen_new = list()
    count = 0
    for m,item in enumerate(resi_index_list_temp):
        if item[1] not in aa_dict.keys():
            if 'X' == seq_antigen[count]:
                seq_antigen_new.append('X')
                resi_index_list_all.append(item[0])
                count += 1
        else:
            if aa_dict[item[1]] == seq_antigen[count]:
                seq_antigen_new.append(aa_dict[item[1]])
                resi_index_list_all.append(item[0])
                count += 1
        if count >= len(seq_antigen):
            break
    seq_antigen_new = ''.join(seq_antigen_new)
    
    if seq_antigen_new != seq_antigen:
        print("Something is wrong with antigen sequence")
    
    #Select pocket indexs
    cmd.select('temp','byres chain {} within {} of chain {}'.format(antigen_chain,dis,Hchain))
    cmd.save('./temp.fasta','temp')
    for fa in SeqIO.parse('./temp.fasta',  "fasta"):
        seq_antigen_pockets.append(antigen_chain+':'+str(fa.seq))
    os.remove('./temp.fasta')
    
    resi_index_list = {'list':[]}
    cmd.iterate('temp',"list.append((resi,resn))",space = resi_index_list)
    resi_index_list_temp = [item for item in list(set(resi_index_list['list']))]
    resi_index_list_temp.sort(key = resi_index_list['list'].index)
    resi_index_list = ','.join([str(resi_index_list_all.index(item[0])) for item in resi_index_list_temp if item[0] in resi_index_list_all])
    # resi_index_list = [resi_index_list_all.index(item[0]) for item in resi_index_list_temp if item[0] in resi_index_list_all]
    # seq_antigen_pockets_new = ''.join([seq_antigen[item] for item in resi_index_list])
    binding_site_antigen = resi_index_list


    #################Step3: Output#################
    
    ##Delete items with very little binding sites
    if (not isinstance(binding_site_nanobody,str)) or (not isinstance(binding_site_antigen,str)):
        cmd.delete('all')
        continue

    if len(binding_site_nanobody.split(','))<10 or len(binding_site_antigen.split(','))<10:
        cmd.delete('all')
        continue
    
    ##Delete items with totally same antigen seq and nanobody seq
    all_seq = seq_Hchain + seq_antigen
    if all_seq not in unique_seq:
        unique_seq.append(all_seq)
    else:
        cmd.delete('all')
        continue
        
    
    
    seq_list.append([pdb,Hchain,seq_Hchain,binding_site_nanobody,antigen_chain,seq_antigen,binding_site_antigen,affinity,affinity_method])
    cmd.delete('all')


column=['pdb','nanobody_chain','seq_nanobody','binding_site_nanobody','antigen_chain','seq_antigen','binding_site_antigen','affinity','affinity_method']
output_dir = './all_binding_site_data_{}A.csv'.format(dis)
output = pd.DataFrame(columns=column,data=seq_list)
output.to_csv(output_dir,index = None)




