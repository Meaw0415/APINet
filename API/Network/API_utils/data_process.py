import numpy as np
import pandas as pd
import pepfeature as pep

## Input: Aptamer Sequence 
## Output: vector with k-mer frequency
def count_kmers(read,k=3):
    """Count kmer occurrences in a given read.

    Parameters
    ----------
    read : string
        A single DNA sequence.
    k : int
        The value of k for which to count kmers.

    Returns
    -------
    counts : dictionary, {'string': int}
        A dictionary of counts keyed by their individual kmers (strings
        of length k).

    Examples
    --------
    >>> count_kmers("GATGAT", 3)
    {'ATG': 1, 'GAT': 2, 'TGA': 1}
    """
    # Start with an empty dictionary
    counts = {}
    counts2={}
    counts3={}
    count_new={'AAA':0,'AAU':0,'AAG':0,'AAC':0,
              'AUA':0,'AUU':0,'AUG':0,'AUC':0,
              'AGA':0,'AGU':0,'AGG':0,'AGC':0,
              'ACA':0,'ACU':0,'ACG':0,'ACC':0,
              'UAA':0,'UAU':0,'UAG':0,'UAC':0,
              'UUA':0,'UUU':0,'UUG':0,'UUC':0,
              'UGA':0,'UGU':0,'UGG':0,'UGC':0,
              'UCA':0,'UCU':0,'UCG':0,'UCC':0,
              'GAA':0,'GAU':0,'GAG':0,'GAC':0,
              'GUA':0,'GUU':0,'GUG':0,'GUC':0,
              'GGA':0,'GGU':0,'GGG':0,'GGC':0,
              'GCA':0,'GCU':0,'GCG':0,'GCC':0,
              'CAA':0,'CAU':0,'CAG':0,'CAC':0,
              'CUA':0,'CUU':0,'CUG':0,'CUC':0,
              'CGA':0,'CGU':0,'CGG':0,'CGC':0,
              'CCA':0,'CCU':0,'CCG':0,'CCC':0,
              }
    # Calculate how many kmers of length k there are
    num_kmers = len(read) - k + 1
    # Loop over the kmer start positions
    for i in range(num_kmers):
        # Slice the string to get the kmer
        kmer = read[i:i+k]
        # Add the kmer to the dictionary if it's not there
        if kmer not in counts:
            counts[kmer] = 0
        # Increment the count for this kmer
        counts[kmer] += 1/num_kmers

    
    for i in counts.keys():
        if i[0]== 'X':
            key1='A'+i[1:]
            key2='C'+i[1:]
            key3='G'+i[1:]
            key4='U'+i[1:]
            counts2[key1]=counts[i]/4
            counts2[key2]=counts[i]/4
            counts2[key3]=counts[i]/4
            counts2[key4]=counts[i]/4
        else:
            counts2[i]=counts[i]

    for i in counts2.keys():
        if i[1]== 'X':
            key1=i[0]+'A'+i[2]
            key2=i[0]+'C'+i[2]
            key3=i[0]+'G'+i[2]
            key4=i[0]+'U'+i[2]
            counts3[key1]=counts2[i]/4
            counts3[key2]=counts2[i]/4
            counts3[key3]=counts2[i]/4
            counts3[key4]=counts2[i]/4
        else:
            counts3[i]=counts2[i]    
    for i in counts3.keys():
        if i[2]== 'X':
            key1=i[0]+i[1]+'A'
            key2=i[0]+i[1]+'C'
            key3=i[0]+i[1]+'G'
            key4=i[0]+i[1]+'U'
            count_new[key1]+=counts3[i]/4
            count_new[key2]+=counts3[i]/4
            count_new[key3]+=counts3[i]/4
            count_new[key4]+=counts3[i]/4
        else:
            count_new[i]=counts3[i]              
    # Return the final counts
    # for key in counts.keys():
    #     counts[key] /= num_kmers

    # count_new["Name"] = name
    
    aptamer_list = []
    for val in count_new.values():
        aptamer_list.append(val)

    return aptamer_list



## Input:  Protein sequence
## Output: Vector with protein feature
def protein_feature(seq):
    protein_list=[]
    
    df = pd.DataFrame(['AKVLVVLLLFAGVDAETHVTGGSAAHAASTFAGLFSPGAKQDIQLINTNGSWHINRTA'],columns=['Info_window_seq'])
    feat = pep.aa_all_feat.calc_df(dataframe=df, aa_column='Info_window_seq', Ncores=1, k=2)
    a = feat.values.tolist()
    
    return a[0][1:]
     
def merge_apt_pro(test_list1,test_list2):
    res_list = [y for x in [test_list1, test_list2] for y in x]
    
    return res_list