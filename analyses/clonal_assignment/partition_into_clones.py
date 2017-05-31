# Code for partitioning Partis output by Daniel Zinder
# Outputs fasta files for different clones, and csv files with their annotations.

# from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
import csv
import sys
import os
from os import makedirs
from tabulate import tabulate
from numpy import histogram

def main(argv):

    #options
    input_fasta_filename = '../../data/Wu_2015_VRC01/Wu_2015_VRC01_heavy_curated.fasta'
    input_partis_csv_filename = '../../results/clonal_assignment/VRC01/VRC01_partition.csv'
    output_filename_stem = 'VRC01_'
    output_fasta_filename_stem = "VRC01_"
    input_annotation_csv_filename = '../../results/clonal_assignment/VRC01/VRC01_annotation.csv'

    remove_letter_list = [] #['.','n','N','?']

    #increase field size for csv
    maxInt = sys.maxsize
    decrement = True

    while decrement:
        # decrease the maxInt value by factor 10 
        # as long as the OverflowError occurs.

        decrement = False
        try:
            csv.field_size_limit(maxInt)
        except OverflowError:
            maxInt = int(maxInt/10)
            decrement = True
        
    #read partitions from csv file
    with open(input_partis_csv_filename) as input_csv_handle:
        partition_reader = csv.reader(input_csv_handle, delimiter=',')
        parsed_csv=[row for row in partition_reader]
        input_csv_handle.close()
        #print(tabulate([row[0:-1] for row in parsed_csv]))
    
    #find best partition    
    partitions=[row[-1] for row in parsed_csv[1:]]
    logprobs=[float(row[0]) for row in parsed_csv[1:]]
    nclusters=[row[1] for row in parsed_csv[1:]]

    best_partition_index=logprobs.index(max(logprobs));
    best_partition=partitions[best_partition_index]
    best_partition_clusters=int(nclusters[best_partition_index])
    print("best partitions with "+ str(best_partition_clusters) +" clusters (logprob="+str(max(logprobs))+"):")
    clone_index=1
    clone_size_list=[]
    for clone in best_partition.split(";"):
        clone_size=len(clone.split(':'))
        clone_size_list.append(clone_size)
        print("clone "+str(clone_index)+" (size "+" "+str(clone_size)+"):\n"+clone+"\n")
    
        clone_index+=1
    
    #parse fasta file
    with open(input_fasta_filename, 'r') as input_fasta_handle:
        input_fasta_dict = SeqIO.to_dict(SeqIO.parse(input_fasta_handle, "fasta"))
        input_fasta_handle.close()
    

    #write clone files
    clone_index=1
    for clone in best_partition.split(";"):
        clone_records=[]
        print("saving clone "+str(clone_index)+" to "+output_filename_stem+str(clone_index)+".fa")
        count=0
        for seqid in clone.split(":"):
            count=count+1
            if count<4:
                print(seqid)
            elif count==5:
                print("....")
            seq_record=input_fasta_dict[seqid]        
            # remove gaps
            remove_gaps=''.join([nuc for nuc in seq_record.seq if nuc not in remove_letter_list])
            seq_record.seq=Seq(remove_gaps)
            clone_records.append(seq_record)
    
        print()
        #os.makedirs(output_fasta_filename_stem,mode=0o777,exist_ok=True) 
        with open('../../results/clonal_assignment/VRC01/'+output_fasta_filename_stem+str(clone_index)+".fasta", 'w') as output_fasta_handle:
            SeqIO.write(clone_records, output_fasta_handle, "fasta")
            output_fasta_handle.close()
    
        clone_index+=1
    
    #read annotation from csv file
    cdr3_length_list=[]
    v_gene_list=[]
    d_gene_list=[]
    j_gene_list=[]
    with open(input_annotation_csv_filename) as input_csv_handle:
        partition_reader = csv.reader(input_csv_handle, delimiter=',')
        parsed_csv=[row for row in partition_reader if len(row[4])>0] # TODO: make skip empty lines 
        input_csv_handle.close()
    
        # for stat reporting (only)
        print("first 20 annotations")
        print(tabulate([row[0:5] for row in parsed_csv[0:20]]))
        cdr3_length_list=[int(row[4]) for row in parsed_csv[1:-1]]
        v_gene_list=[row[1] for row in parsed_csv[1:-1]]
        d_gene_list=[row[2] for row in parsed_csv[1:-1]]
        j_gene_list=[row[3] for row in parsed_csv[1:-1]]
    
    v_dict={row[0]:row[1] for row in parsed_csv[1:-1]}
    d_dict={row[0]:row[2] for row in parsed_csv[1:-1]}
    j_dict={row[0]:row[3] for row in parsed_csv[1:-1]}
    cdr3_length_dict={row[0]:row[4] for row in parsed_csv[1:-1]}

    # write clone annotation into seperate files   
    clone_index=1
    for clone in best_partition.split(";"):
        clone_records=[]
        print("saving clone "+str(clone_index)+" to "+output_filename_stem+str(clone_index)+".csv")
    
    
    
        with open('../../results/clonal_assignment/VRC01/'+output_fasta_filename_stem+str(clone_index)+".csv", 'w') as output_csv_handle:
    
            vdjwriter = csv.writer(output_csv_handle, delimiter=',')
        
            vdjwriter.writerow(['id','cdr3len','v','d','j'])
            {vdjwriter.writerow([seqid,cdr3_length_dict[seqid],v_dict[seqid],d_dict[seqid],j_dict[seqid]]) for seqid in clone.split(":") if seqid in v_dict.keys()}  
            output_csv_handle.close()
        
  
        clone_index+=1

if(__name__ == "__main__"):
    status = main(sys.argv)
    sys.exit(status)



