## Alignments ##

We used [MACSE](http://bioweb.supagro.inra.fr/macse/) to align observed sequences along with a "germline" sequence consisting of the concatenated V and J genes. For the light chains, we used only V genes (V + J seems to include frameshifts from the germline sometimes).

We refer to each of the seven heavy and light chain datasets as a "clone": CH103, CH103L, VRC26, VRC26L, VRC01-13, VRC01-01, VRC01-19. The "L" in CH103L and VRC26L indicates that those clones are the light chains of lineages CH103 and VRC26.

The steps used to process and align original fasta files (and fasta files obtained after partitioning the VRC01 datasets) are detailed below. However, the final alignment files are provided in the ```results``` directory of this repository as ```<CLONE>_final_alignment.nex```. 

### Steps ###
0 - If the user has the ```MACSE.jar``` file in his or her home folder (in a Unix system), executing ```MACSE_alignment.sh``` will execute steps 1-2 for all datasets.

1 -  In the ```analyses/alignments/``` directory, execute the python scripts ```prepare_<CLONE>_for_MACSE.py``` to add a concatenated germline sequence to the fasta files with the observed sequences (```<CLONE>_plus_GERMLINE.fasta```). 

2 - Run MACSE: 

```java -jar path/to/MACSE/MACSE.jar -prog alignSequences -seq ../../results/alignments/<CLONE>_plus_GERMLINE.fasta -../../results/alignments/<CLONE>_MACSE_NT.fasta -../../results/alignments/<CLONE>_MACSE_AA.fasta```

2 - Manually edit nucleotide alignments (```<CLONE>_MACSE_NT.fasta```) and save in nexus format as ```<CLONE>_MACSE_NT_EDITED.nex```. For heavy-chain alignments, make sure there are no gaps within the V and within the J portion of the "GERMLINE" sequence (the ones that are concatenated to produce the alignment). There should be a large gap between them (corresponding to the D gene from the full VDJ sequences from the clones). We removed mostly-gap regions at the edges of the alignment.

3 - Run ```add_dates_<CLONE>.py``` scripts to execute python scripts that will add sequence sampling times from the original fasta files to the edited alignment files and output the final alignment files (```<CLONE>_final_alignment.nex```)

### Deleted sequences ###

During the manual edition step, we deleted the following sequences:

**Alignment of VRC26 light chain**: removed sequences that contained frameshifts in the MACSE alignment; those sequences are probably incomplete reads with a few "extra" or missing nucleotides at the beginning of the sequence. We removed the following sequences (9 out of 472 sequences, approx. 2%) from the MACSE alignment and did not the remaining sequences again (identified by their GenBank accession numbers):

KJ134566, KJ134527 ,KJ134602, KJ134521, KJ134482, KJ134504, KJ134594, KJ134846, KJ134750

**Alignment of VRC01-19** - removed one sequence with frameshift indels in MACSE alignment (KP841240).
