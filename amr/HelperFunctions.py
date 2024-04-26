import requests as r
from Bio import SeqIO
from io import StringIO
import pandas as pd
import re


'''takes in a protein id, cID, and returns the complete amino acid sequence''' 
def get_protein_seq(cID):
    baseUrl="http://www.uniprot.org/uniprot/"
    currentUrl=baseUrl+cID+".fasta"
    response = r.post(currentUrl)
    cData=''.join(response.text)
    
    Seq=StringIO(cData)
    pSeq=list(SeqIO.parse(Seq,'fasta'))

    return str(pSeq[0].seq)

'''Takes in an amino acid, AA, and peptide sequence, PEPTIDE, which are both strings. 
Finds the indices within the peptide where there is a modifification. PEPTIDE string must be of one
of the following formats:

    AADTIGYPVM[649.3660]IR
    AADTIGYPVM(49.3660)IR

The function returns the zero indexed position of the peptide if there is a single amino acid of type 
AA modified. Returns -1 if there are multiple amino acids of type AA modified or no amino acids modified.
   ''' 
def find_position_peptide_sequence(aa, peptide):
    


    return 

"""This function takes in a modified peptide sequence, modified_peptide_seq, an amino acid, amino_acid, and 
a complete protein sequence, complete_protein_sequence. modified_peptide_seq is of the form, IAM[649.3660]QTLDMGR, where
M is the modified amino acid and [649.3660] is the modification. This function finds the index in the complete_protein_sequence
in which the modification appears."""
def find_the_index_of_the_modification(modified_peptide_seq, amino_acid, complete_protein_sequence):
    # Extract the modification and its index in the modified peptide sequence
    modification = re.search(fr'{amino_acid}\[(.*?)\]', modified_peptide_seq)
    print(modification.groups())

    # mod_index = modified_peptide_seq.index('[' + modification + ']')

    # # Find the index of the amino acid in the complete protein sequence
    # amino_acid_index = complete_protein_sequence.find(amino_acid)

    # # Find the offset caused by any modifications before the amino acid
    # offset = mod_index - modified_peptide_seq.index(amino_acid)

    # # Calculate the index in the complete protein sequence
    # index_in_protein = amino_acid_index + offset

    # return index_in_protein  