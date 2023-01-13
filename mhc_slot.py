# MHC_slot_dataset
# Shalabh Shukla

# The purpose is to construct different molecular configurations of binding slots for alleles of interest using RDkit, returns smile 
# data as alternate x value instead of AAsq. The three working models that will be tested include: 'Full_Binding_Slot', 'AASeq_2Helices',
#'AASeq_2Helices_Cyclic', can include residue models / other models if needed. Match y values to data using Allele names (allele name     
# nomencalature if finicky see allele entry note below). 1 csv file  required for code to work if no pickles present, one sequnce database 
# file for allknown MHC sequences (update for new species / isoforms.

""" Hyperparameters for narrowing search and alignment (tune for use cases) """

#BINDING_DOMAIN_REF = [
#   "SNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQ"     # Sars Cov-2 binding domain of spike protein
#    ]

BINDING_DOMAIN_REF = [
    "PWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRG",                                       # Refrence Sequnece helix 1 for B:35:01 isoform AASeq[73:107] 
    "AQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKET"                                   # Refrence Sequnece helix 2 for B:35:01 isoform AASeq[163:202]
    ]

#BINDING_DOMAIN_REF = [
#   "GRFASFEAQGALANIAVDKANLEIMTKRSNYTPITNVPPEVTVLTNSPVELR",                     # Refrence Sequnece helix 1 for HLA-DRA1 
#   "TELGRPDAEYWNSQKDLLEQRRAAVDTYCRHNYGVG"                                      # Refrence Sequnece helix 2 for HLA-DRB1
#   ]

THRESHOLD_COVERAGE = 0.30                                                       # percent alignment for a succesful match (0-1 scale)

G_LINKER_LENGTH = 4                                                             # varies depending on size of binding slot of the protein

""" Needed modules for mhc_slot_dataset """

import os
import csv
import pickle
import pandas as pd 
from rdkit import Chem
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62                                     # matrix optimizes AA reads in search algorithm

def mhc_slot_dataset(
    allele,
    binding_model,
    return_type,
    sequences_database = 'MHC_Full_AAseq_DB.csv',                               # sequences_database = 'MHC_Full_AAseq_DB.csv' repository of all sequences and associated allele names
    directory = '../data/alleles/'):
    """ mhc_slot_dataset

    inputs:
    allele - protein variant name included in sequences_database
    
    binding_model - various models of protein regions that participate in
    binding such as: full_complex which is simply the full protein sequence,
    binding slot which is the full region start to end that is participating in
    binding, linear which is disparate binding regions linkead linearly, cyclic
    which is disparate binding regions linkead cyclically

    return_type - data generated as either stringa (amino acid sequence) or
    smiles outputs

    sequences_database - csv file which contains allele names in one column and
    accompanying amino acid sequence in a second column

    directory - location of folder that contains the sequences_database csv file 
    """
   
    if not (binding_model == 'full_complex') | \
           (binding_model == 'binding_slot') | \
           (binding_model == 'linear') | \
           (binding_model == 'cyclic'):
        
        raise ValueError('binding_model must be: full_complex, \
                          binding_slot, helices, or cyclic')     
        
    if not (return_type == 'strings') | (return_type == 'smiles'):
        
        raise ValueError('return_type must be: strings or smiles')  
                  
    
    def full_seq_2_domain_seq(seq_full,seq_domain):
        """ Helix match function 
    
        inputs:
        sequence read
        refrence sequences
        """
        
        aligned_region = 0; domain_adjustment = 0                               # helix allignemnt, pairwise needs 5 args returns 5 elements in list
        
        aligned_domain = pairwise2.align.localds(
            seq_full, seq_domain, blosum62, -10, -0.5)[0][-3:]
        
        domain_diferences = (len(seq_domain)
            - (aligned_domain[2] - aligned_domain[1]))                          # alignment starts and ends on matches only
        
        if domain_diferences > 0:                                               # allignment not contigous to refrence helix if diffrence greater than 0
            
            for domain_position in range(len(seq_domain) - 1):
                
                if seq_full[aligned_domain[1]] == \
                    seq_domain[domain_position] or \
                    seq_full[aligned_domain[1] + 1] == \
                    seq_domain[domain_position + 1]:
                        
                    aligned_region = seq_full[
                        aligned_domain[1] - domain_position:
                        aligned_domain[2] + domain_diferences - domain_position
                        ]                                                       # allignment adjustments for mismatch
                    
                    domain_adjustment = (aligned_domain[2] 
                        + domain_diferences 
                        - domain_position)
                    
                    return [aligned_region, aligned_domain, domain_adjustment] 
                     
        else:
            
            aligned_region = seq_full[aligned_domain[1]:aligned_domain[2]]
            
            domain_adjustment = aligned_domain[2]                                # perfect allignment case / insertion case
        
        return [aligned_region, aligned_domain, domain_adjustment]
    
    def domain_seq_2_string_type(seq_full, seq_aligned):
        """ Return type data generation function
    
        inputs:
        sequence read 
        aligned sequence with refrence domain 
        """
        
        AAseq_full = seq_full; AASeq_slot = '' 
        AASeq_domain_lin = '' ; AASeq_domain_cyc = ''
        AASeq_slot = seq_full[:seq_aligned[-1][2]]                              # full Sequence of binding slot, multi-complex case returns combined sequence 
        
        for refrences in range(len(BINDING_DOMAIN_REF)): 
            AASeq_domain_cyc = (AASeq_domain_cyc 
                + seq_aligned[refrences][0] 
                + ("G" * G_LINKER_LENGTH))                                      # cyclic linkage of binding domains
        
        AASeq_domain_lin = AASeq_domain_cyc[:-G_LINKER_LENGTH]                  # linear linkage of binding domains
        
        return [AAseq_full, AASeq_slot, AASeq_domain_lin, AASeq_domain_cyc]     # string sequences of return type peptides     
    
    def string_type_2_smile_type(binding_model_seqs):
        """ Smiles data generation for all return types function 

        inputs:
        sequences of specified binding model
        """
        
        binding_model_smiles = binding_model_seqs
        
        for types in range(len(binding_model_seqs) - 1):
            
            binding_model_seqs[types] = Chem.MolToSmiles(
                Chem.MolFromHELM("PEPTIDE1{" 
                    + ('.'.join(binding_model_seqs[types])) 
                    +"}$$$$V2.0" ))                                             # SMILES output  
                                                                                                      
        binding_model_smiles[-1] = Chem.MolToSmiles(
            Chem.MolFromHELM("PEPTIDE1{" 
                + '.'.join(binding_model_seqs[-1]) 
                + "}$PEPTIDE1,PEPTIDE1,1:R1-" 
                +  str(len(binding_model_seqs[-1])) 
                + ":R2$$$V2.0"))                                                # cyclic HELM
        
        return binding_model_smiles 
    
    def alignment_graph(seq_full, BINDING_DOMAIN_REF): 
        """ Show alingment on failure """
    
        for refrences in range(len(BINDING_DOMAIN_REF)):
            print(BINDING_DOMAIN_REF[refrences] + '\n')
            
            matches = pairwise2.align.localds(
                seq_full, BINDING_DOMAIN_REF[refrences], blosum62, -10, -0.5)   # start sequence for HLA of interest
            
            for alignment in matches:
                print(pairwise2.format_alignment(*matches[0]))
    
    def pickle_return_type(model_data):
        """ Pickle requested return_type data """
        
        protein_return_type_df = pd.DataFrame(
                [{
                    'full_complex': model_data[0],
                    'binding_slot': model_data[1],
                    'linear': model_data[2],
                    'cyclic': model_data[3]
                }],
                index = [allele_serotype],
                columns = ['full_complex', 'binding_slot', 'linear', 'cyclic']) # store binding_models in dataframe
        
        pd.to_pickle(
            protein_return_type_df, 
            os.path.join(directory, 
                (allele_name + '_' + return_type + '.pickle')))                 # generate pickle 
        
        return protein_return_type_df
                
    """ Check if allele is already pickled in directory 
    if not generate new dataset or allele pickle from new dataframe 
    """                 
    
    try:               
                                                                                # multi-complex case 
        allele_serotype = ''.join(allele)  
    
    except:          
                                                                                # single allele case
        allele_serotype = str(allele)
    
    allele_name = allele_serotype
        
    for char in ['*', ':']:  
                                                                                # for some reason they decided to use math operators for their allele naming convention
        allele_name = allele_name.replace(char, "")                             # removing chars still result in unique identifiers for alleles of interest

    if os.path.isfile(
        os.path.join(directory,(allele_name + '_' + return_type + '.pickle'))): # check for allele pickle
        
        print('\nLoaded ' 
            + binding_model 
            + ' data for ' 
            + allele_serotype 
            + ' from ' 
            + (allele_name + '_' + return_type + '.pickle \n'))
        
        protein_return_type_df = pickle.load(
            open(
                os.path.join(
                    directory, 
                    (allele_name + '_' + return_type + '.pickle')),
                "rb"))                                                          # load from allele pickle
        
        return (protein_return_type_df.loc[allele_serotype, binding_model])     # returns values from pickle
        
    """ Generating dataset for allele of interest from sequences csv """ 
    
    try:
        
        file = open(os.path.join(directory,sequences_database)); file.close() 
        
    except FileNotFoundError:
        
        print(sequences_database + ' does not exist or not in ' + directory)
        raise
        
    print('\nMaking ' + return_type +  ' pickle for ' + allele_serotype + '\n')
    aligned_domain = []; threshold_score = []                                   # intizing allignments and score thresholds
    
    for refrences in range(len(BINDING_DOMAIN_REF)):
        
        threshold_score.append((
            pairwise2.align.localds(
                BINDING_DOMAIN_REF[refrences],
                BINDING_DOMAIN_REF[refrences],
                blosum62, -10, -0.5))[0][2]
            * THRESHOLD_COVERAGE)
    
    """ Sequences search and stitching together disparate domains 
    to form full binding slot sequence across multpile complexes / regions """   
    
    with open (os.path.join(directory, sequences_database), 'r') as seq_csv:
            
        seq_full = []; subunits = 0 
        sequence_reads = csv.reader(seq_csv) 
        
        for sequence_read in sequence_reads:
            
            if subunits < len(allele) and \
                allele[subunits] == sequence_read[0][0:len(allele[subunits])]:  # if multiple subunits specified stitches them together to form single binding slot sequence
                
                seq_full.append(sequence_read[1]); subunits += 1 
            
        for refrences in range(len(BINDING_DOMAIN_REF)):
            
            aligned_domain.append(
                full_seq_2_domain_seq(
                    ''.join(seq_full),
                    BINDING_DOMAIN_REF[refrences]))                             # search for aligned sequences to refrence domains
            
            if (aligned_domain[refrences])[1][0] < threshold_score[refrences]:  # aligned region fails minimum threshold score
                
                alignment_graph(''.join(seq_full), BINDING_DOMAIN_REF)          # return alignment graph to show mismatches
                raise NameError('Input allele did not allign with refernce')
                  
        binding_protein_str = domain_seq_2_string_type(
            ''.join(seq_full), aligned_domain)
        
        if return_type == 'strings':
            
            protein_return_type_df = pickle_return_type(binding_protein_str)
        
        if return_type == 'smiles':
            
            binding_protein_smiles = string_type_2_smile_type(
                binding_protein_str)
            
            protein_return_type_df = pickle_return_type(binding_protein_smiles)
                    
        return (protein_return_type_df.loc[allele_serotype, binding_model])
    
    raise NameError('input allele ' 
        + allele_serotype + ' not in ' 
        + sequences_database)
