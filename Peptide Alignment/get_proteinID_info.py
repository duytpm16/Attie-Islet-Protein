
# coding: utf-8

# In[2]:


import csv
import pandas
import pickle
import os


# In[6]:


blast_results = pandas.read_csv('blast_results.txt', sep = '\t', low_memory = False)


# In[7]:


pkl_file = open('mm_peptide_info.pkl', 'rb')
peptide_info = pickle.load(pkl_file)
pkl_file.close()


# In[8]:


blast_results


# In[9]:


peptide_info


# In[34]:


subject_gene_id       = []
subject_gene_symbol   = []
subject_transcript_id = []
subject_chromosome    = []




for i in range(0, len(blast_results)):
    if not pandas.isnull(blast_results['subjects'][i]):
        subject = blast_results['subjects'][i].split(',')
        
        subject_gene_id.append(','.join(peptide_info[s]['gene_id'] for s in subject))
        subject_gene_symbol.append(','.join(peptide_info[s]['gene_symbol'] for s in subject))
        subject_transcript_id.append(','.join(peptide_info[s]['transcript_id'] for s in subject))
        subject_chromosome.append(','.join(peptide_info[s]['chromosome'].split(':')[1] for s in subject))
        
    else:
        
        subject_gene_id.append('NA')
        subject_gene_symbol.append('NA')
        subject_transcript_id.append('NA')
        subject_chromosome.append('NA')


# In[35]:


subject_gene_id       = pandas.Series(subject_gene_id)
subject_gene_symbol   = pandas.Series(subject_gene_symbol)
subject_transcript_id = pandas.Series(subject_transcript_id)
subject_chromosome    = pandas.Series(subject_chromosome)


# In[36]:


blast_results['gene_id'] = subject_gene_id
blast_results['gene_symbol'] = subject_gene_symbol
blast_results['transcript_id'] = subject_transcript_id
blast_results['chromosome'] = subject_chromosome 


# In[37]:


blast_results


# In[39]:


blast_results.to_csv('blast_results.txt', sep = '\t', index = False, header = True)

