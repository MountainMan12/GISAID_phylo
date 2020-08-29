import pandas as pd
import Bio
from Bio import SeqIO
from functools import reduce

#Reading the clinical data
clinical_data2 = pd.read_csv('gisaid_hcov-19_table.tsv', delimiter ='\t')
#adding undersccore to column names
clinical_data2.columns = [c.replace(' ','_') for c in clinical_data2.columns]
#Drop NAs in patient age and status columns
clinical_data_filt2 = clinical_data2.dropna(how='any', subset =['Patient_status'])
clinical_data_filt3 = clinical_data_filt2.dropna(how='any', subset =['Patient_age'])
#Removing samples with unknown patient status
missing1 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'unknown ']
missing2 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'Not known']
missing3 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'unkown']
missing4 = clinical_data_filt3[clinical_data_filt3.Patient_status == ' unknown']
missing5 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'unknow']
missing6 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'Unknow']
missing7 = clinical_data_filt3[clinical_data_filt3.Patient_status == '\ufeffunknown']
missing8 = clinical_data_filt3[clinical_data_filt3.Patient_status == '-']
missing9 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'uncknown']
missing10 = clinical_data_filt3[clinical_data_filt3.Patient_status == 'Unkown']

#Removing samples for non human host
missing11 = clinical_data_filt3[clinical_data_filt3.Host == 'Environment']
missing12 = clinical_data_filt3[clinical_data_filt3.Host == 'Panthera tigris jacksoni']

#Removing sample with bad sequence
missing13 = clinical_data_filt3[clinical_data_filt3.Accession_ID == 'EPI_ISL_494759']

#Removing samples with unknown patient age
missing14 = clinical_data_filt3[clinical_data_filt3.Patient_age == 'unknown']
missing15 = clinical_data_filt3[clinical_data_filt3.Patient_age == 'Unknown']

#Removing samples with bad age 
missing16 = clinical_data_filt3[clinical_data_filt3.Patient_age == 'unkown']

#Removing all samples corresponding to missing values
ind = []
for i in range(1,17):
    ind.append(list(globals()["missing" + str(i)].index)) #Appending all indexes in a list
#Conversion to a 1-D list
single_list = reduce(lambda x,y: x+y, ind) 

single_set = set(single_list) #Unique indices
clinical_data_filt4 = clinical_data_filt3.drop(single_set)

#Writing to a new file
clinical_data_filt4.to_csv('gisaid_hcov-clinical_filt.tsv', sep ='\t', index = False)

new_ids =list(clinical_data_filt4["Accession_ID"])

#Writing the filtered fasta sequences using the clinical accession IDs
new_file = open('gisaid_hcov-19_filt.fasta','w')
j = 0
for records in SeqIO.parse('gisaid_hcov-19_2020_08_17_14.fasta','fasta'):
    for i in new_ids:
        if i in records.description:
            j += 1
            new_file.write('>'+ i + ' | ' + records.description.split('/')[1] + str(j))
            new_file.write('\n')
            new_file.write(str(records.seq))
            new_file.write('\n')

new_file.close()