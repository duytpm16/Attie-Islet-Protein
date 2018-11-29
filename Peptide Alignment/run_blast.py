### Load Libraries
from Bio import SeqIO
import subprocess
import csv
import pandas
import os







### Read in MaxQuant Output File
data = pandas.read_csv('peptides_DO_islets.txt', sep = '\t', low_memory = False)






# ### Dimension of MaxQuant File
n_rows = data.shape[0]
data.shape



### Remove ROWS in MaxQuant File where SEQUENCES are not Peptides
amino_acids = ['G','A','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
rows_to_drop = []


for i in range(0, n_rows):
    if data['Sequence'][i][0] not in amino_acids:
        rows_to_drop.append(i)

        
data = data.drop(rows_to_drop)
data.index = range(len(data))
data.shape



### For each Peptide Sequence:
#      1.) Create temporary fasta file
#      2.) Do blast alignment
#      3.) Find subjects that have 100% identity alignment length match to the query peptide
peptide_sequence = []
identity_percent = []
blast_subject    = []
align_length     = []
subject_start    = []
subject_end      = []

for i in range(0,len(data['Sequence'])):
    
    # 1.) Create the temporary fasta file
    with open('output_fasta.fa', 'w') as out_file:
        header = '>' + 'row' + str(i) + '\n'
        out_file.write(header)
        out_file.write(data['Sequence'][i] + '\n')
    
    
    
    
    
    # 2.) Blast alignment begins
    subprocess.call(['/usr/local/ncbi/blast/bin/blastp', '-db', 'Mm_pep_blast', '-query', 'output_fasta.fa' ,'-comp_based_stats','0', '-outfmt', '10', '-out','test'])
    
    
    
    
    
    
    # 3.) If there are subjects that the query aligned to, do..
    if os.stat("test").st_size != 0:
        
        
        # Read in file. Give column names since blast alignment output doesn't give
        #    Note*: change -outfmt to 7 in 2.) to see header names
        blast_output = pandas.read_csv('test', sep = ',', header = None)
        blast_output.columns = ['query', 'subject', 'identity', 'alignment_length', 'mismatches', 'gap_opens', 'q.start', 'q.end', 's.start','s.end', 'evalue', 'bit_score']
    
    
    
    
    
    
    
    
        # Find alignments that match requirement
        identity_100 = []
        for j in range(len(blast_output)):
            if (blast_output['identity'][j] == 100.0) and (blast_output['alignment_length'][j] == len(data['Sequence'][i])):
                identity_100.append(j)
    
    

    
    
    
    
    
    
        # Subset blast_output to those that match 3.)
        blast_output = blast_output.iloc[identity_100,]
    
     
    
    
        
    
        
        
        
        
        #If there are still subjects after filter, then keep, else insert NA
        if(len(blast_output) != 0):
            peptide_sequence.append(data['Sequence'][i])
            blast_subject.append(','.join(str(s) for s in blast_output['subject']))
            identity_percent.append(','.join(str(s) for s in blast_output['identity']))
            align_length.append(','.join(str(s) for s in blast_output['alignment_length']))
            subject_start.append(','.join(str(s) for s in blast_output['s.start']))
            subject_end.append(','.join(str(s) for s in blast_output['s.end']))
        else:
            peptide_sequence.append(data['Sequence'][i])
            blast_subject.append('NA')
            identity_percent.append('NA')
            align_length.append('NA')
            subject_start.append('NA')
            subject_end.append('NA')
    
    
    
    
    
    else:    
        
        # If no alignments append NA
        peptide_sequence.append(data['Sequence'][i])
        blast_subject.append('NA')
        identity_percent.append('NA')
        align_length.append('NA')
        subject_start.append('NA')
        subject_end.append('NA')
    



    
    
    
    
    
    
### Save blast results to data.frame
blast_result = pandas.DataFrame({'peptide_sequence': peptide_sequence,
                                 'subjects':         blast_subject,
                                 'identity_percent': identity_percent,
                                 'alignment_length': align_length,
                                 'subject_start':    subject_start,
                                 'subject_end':      subject_end})
   



    
    
    
    
    
    
    
### Save dataframe as .txt
blast_result.to_csv('blast_results.txt', sep = '\t', index = False, header = True)

