#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Imports
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import MeltingTemp as mt
import pandas as pd


# In[ ]:


def fasta_input(probe_complete, gene_name):
    '''Converts the probe output into a csv convenient for FASTA inputs'''
    # Put the sequences in a list
    bip = [i for i in probe_complete.iloc[:,1]]
    
    # Match the sequences to their names
    for i in range(len(probe_complete)):
        bip.insert(2*i, '>' + probe_complete.iloc[i,0])
        
    # Convert list to a pandas dataframe
    bip_pd = pd.DataFrame({'Probes':bip})
    
    # Name the file
    name = gene_name + '_fasta-input.csv'
    
    # Export
    bip_pd.to_csv(name)
    
    return bip_pd


# In[ ]:


def rap_probes(file, gene_name, probe_length = 50, coverage  = 1, adaptor = 'CAAGTCA', nmol = '25nm'):
    '''Takes a FASTA imported as a Pandas DataFrame and generates 
    RNA affinity purification probes based on given parameters'''
    # Necessary imports:
    # from Bio.Seq import Seq
    # import pandas as pd
    # from Bio import SeqIO
    # from Bio.SeqUtils import MeltingTemp as mt
    # fasta_input function defined above
    
    # Import the FASTA as a Pandas Dataframe
    ls = SeqIO.read(file, 'fasta').seq
    ls = ls.reverse_complement()
    
    # Generate start indices for each probe
    ind2 = [int(i*(probe_length/coverage)) for i in range(0,int(len(ls)/probe_length*coverage))]
    
    # Add a last probe if more than a quarter probe remains
    if ind2[-1] < int(len(ls) - probe_length/4):
        ind2.append(len(ls)-probe_length)
    
    # Generate the probes
    probes = [ls[i:i+probe_length] for i in ind2]
    
    # Add adaptor to 5' end of probes
    probes = [str(adaptor + i) for i in probes]
    
    # Double check probe lengths
    lengths = [len(i) for i in probes]
    
    # Calculate Tms for probes based on hybridization buffer
    melts = [mt.Tm_NN(i, Na = 500, Tris = 10) for i in probes]
    melts = [round(i,1) for i in melts]
    
    # Make appropriate labels for each probe
    # Saving name as a second variable makes later design steps easier
    name = gene_name + '_' + str(probe_length) + 'mer_'
    labels = [name + str(i) for i in range(1, len(probes)+ 1)]
    
    # Add nanomoles to order from IDT
    idt_bulk = [nmol for i in range(len(probes))]
    
    # Put everything in a DataFrame
    probe_complete = pd.DataFrame({
        'Name':labels,
        'Sequence':probes,
        'Length (bp)': lengths,
        'Tm (°C)':melts,
        'Nanomoles':idt_bulk
    })
    
    # Name your file
    filename = gene_name + '_' + str(probe_length) + 'mer' + '_Probes.csv'
    
    probe_complete.to_csv(filename)
    
    
    # Create a csv that is compatible with programs that require fasta inputs
    fasta_name = gene_name + '_' + str(probe_length) + 'mer' + '_Probes'
    
    fasta_input(probe_complete, fasta_name)
    
    # Export the following:
    # probe_complete: pandas array with the full set of probes
    # filename: name of the file; will be used to update csv after filtering
    # name: header for each probe; useful for indexing probe_complete
    return probe_complete, filename, name


# In[ ]:


def probe_filter(bad_probes, probe_complete, name, filter_type = 'blat'):
    '''Takes a list of probes with undesirable features (>25bp blat match or 
    repeat elements) and filters them from the full list'''
    bad_probes = [name + str(i) for i in bad_probes]
    
    for i in bad_probes:
        probe_complete = probe_complete[probe_complete.Name != i]
    
    filtername = name + filter_type + '-filtered'
    
    fasta_input(probe_complete, filtername)
    
    return probe_complete

# In[ ]:
def reindex(probe_complete, name):
    for i in range(len(probe_complete)):
        probe_complete.iloc[i,0] = name + str(i+1)
    
    return probe_complete
