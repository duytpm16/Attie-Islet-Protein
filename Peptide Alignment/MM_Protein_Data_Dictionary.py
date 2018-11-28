
# coding: utf-8

# ### Load Libraries

# In[20]:


from Bio import SeqIO
import pickle


# <br>
# <br>
# <br>
# <br>
# ### Create data dictionary for each protein in the Mus Musculus Protein Fasta File

# In[25]:


peptide_info = {}

for seq_record in SeqIO.parse('Mus_musculus.GRCm38.pep.all.fa','fasta'):

    temp_description = seq_record.description.split(' ')

    if len(temp_description) > 8:
        peptide_info[temp_description[0]] = {'chromosome'   : temp_description[2].replace('chromosome:',''),
                                             'description'  : temp_description[8].replace('description:',''),
                                             'gene_biotype' : temp_description[5].replace('gene_biotype:',''),
                                             'gene_id'      : temp_description[3].replace('gene:',''),
                                             'gene_symbol'  : temp_description[7].replace('gene_symbol:',''),
                                             'sequence'     : str(seq_record.seq),
                                             'transcript_id': temp_description[4].replace('transcript:','')}
    else:
        peptide_info[temp_description[0]] = {'chromosome'   : temp_description[2].replace('chromosome:',''),
                                             'gene_biotype' : temp_description[5].replace('gene_biotype:',''),
                                             'gene_id'      : temp_description[3].replace('gene:',''),
                                             'gene_symbol'  : temp_description[7].replace('gene_symbol:',''),
                                             'sequence'     : str(seq_record.seq),
                                             'transcript_id': temp_description[4].replace('transcript:','')}
        
        


# <br>
# <br>
# <br>
# <br>
# ### Save the data dictionary

# In[26]:


output = open('mm_peptide_info.pkl', 'wb')

pickle.dump(peptide_info, output)

output.close()

