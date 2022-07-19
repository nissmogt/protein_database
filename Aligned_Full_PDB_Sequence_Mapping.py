#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Author: Hana Jaafari
#Date: 04/27/2021
#Purpose: The purpose of this script is to map aligned, pdb, and full versions of a sequence to each other.
import urllib
import pickle
import numpy as np
from Bio.PDB import *


# In[ ]:


def aligned_pdb_sequence_mapping(protein_identifier,pairwise_aligned_full_pdb_sequence,
                                 pairwise_aligned_pdb_sequence,aligned_protein_sequence,
                                 filtered_aligned_protein_sequence,completed_pdb_directory=None):
    parser = PDBParser()
    if completed_pdb_directory==None:
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{protein_identifier[:4]}.pdb",
                                   f"./{protein_identifier[:4]}.pdb")
        pdb_structure_file=f"./{protein_identifier[:4]}.pdb"
    else:
        pdb_structure_file=f"{completed_pdb_directory}/{protein_identifier}.pdb"
    structure = parser.get_structure(protein_identifier,pdb_structure_file)
    chain=protein_identifier.split("_")[1]
    ###
    complete_to_original_dict={};counter=0
    for i,x in enumerate(pairwise_aligned_full_pdb_sequence):
        if pairwise_aligned_pdb_sequence[i]=="-":
            complete_to_original_dict[i]=np.nan
        else:
            complete_to_original_dict[i]=counter
            counter+=1
    print("Mapping from complete pdb sequence to cleaned sequence is:"); print(complete_to_original_dict)
    pickle.dump(complete_to_original_dict, open(f"{protein_identifier}_complete_to_cleaned_sequence_map.pickle", 'wb'))
    ###
    residue_index_dict={};counter=0
    for i in structure[0][chain].get_residues():
        if i.get_full_id()[3][0].startswith("H_"):
            residue_index_dict[i.get_full_id()[3][1]]=np.nan
        else:
            residue_index_dict[i.get_full_id()[3][1]]=counter
            counter+=1
    print("Mapping from full pdb sequence to cleaned sequence is:"); print(residue_index_dict)
    pickle.dump(residue_index_dict, open(f"{protein_identifier}_full_to_cleaned_sequence_map.pickle", 'wb'))
    ###
    full_to_aligned_index_dict={}; counter=0
    for i,x in enumerate(aligned_protein_sequence):
        if x != "-" and x==x.upper():
            full_to_aligned_index_dict[counter]=i
        if x!="-":
            counter+=1
    print("Mapping from cleaned sequence to aligned sequence is:"); print(full_to_aligned_index_dict)
    pickle.dump(full_to_aligned_index_dict, open(f"{protein_identifier}_cleaned_to_aligned_sequence_map.pickle", 'wb'))
    ###
    dash_indices=[i for i,x in enumerate(filtered_aligned_protein_sequence) if x!="-"]
    print(len(dash_indices))
    counter=0
    for entry in full_to_aligned_index_dict:
        full_to_aligned_index_dict[entry]=dash_indices[counter]
        counter+=1

    print("Mapping from cleaned sequence to filtered aligned sequence is:"); print(full_to_aligned_index_dict)
    pickle.dump(full_to_aligned_index_dict, open(f"{protein_identifier}_cleaned_to_filtered_aligned_sequence_map.pickle", 'wb'))
    ###
    pdb_to_aligned_index_dict={}
    for entry in residue_index_dict:
        if residue_index_dict[entry] in full_to_aligned_index_dict:
            correct_index=residue_index_dict[entry]
            filtered_correct_index=full_to_aligned_index_dict[correct_index]
            pdb_to_aligned_index_dict[entry]=filtered_correct_index
            residue=filtered_aligned_protein_sequence[filtered_correct_index]
        else:
            filtered_correct_index=np.nan
            pdb_to_aligned_index_dict[entry]=filtered_correct_index
            residue="-"
        with open(f"{protein_identifier}_Filtered_Aligned_to_PDB_Sequence_Mapping.txt","a") as f:
            f.write("{} {} {} \n".format(filtered_correct_index,entry,residue))
    print("Mapping from filtered aligned sequence to full pdb sequence is:"); print(pdb_to_aligned_index_dict)
    pickle.dump(pdb_to_aligned_index_dict, open(f"{protein_identifier}_filtered_aligned_to_full_sequence_map.pickle", 'wb')) 
    
    return residue_index_dict,full_to_aligned_index_dict

