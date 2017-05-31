#ALIGNMENTS#

We used MACSE to align observed sequences along with a "germline" sequence consisting of the concatenated V and J genes (as identified by authors in the corresponding papers for each lineage). For the light chains, we used only V genes (V + J seems to include  frameshifts from the germline sometimes)

```*_MACSE_NT.fasta``` files are the raw output files from MACSE.

```*_MACSE_NT_EDITED.nex``` files are the manually adjusted alignment files without sampling time annotation.

```*_final_alignment.nex``` files are the final alignment files with annotated sampling times, ready to be used as input by BEAST.

##Steps##

1 - Run MACSE\_alignment.sh. This will run the python scripts (prepare\_<CLONE>\_for\_MACSE.py) that add a concatenated germline sequence to the fasta files with the observed sequences. Temporary files containing the germline sequence (<CLONE>\_plus\_GERMLINE.fasta>) are then passed on to MACSE to be aligned. MACSE in turn outputs <CLONE>\_MACSE\_<AA/NT>.fasta files.

2 - Manually edit nucleotide alignments in SeaView and save in nexus format as <CLONE>\_MACSE\_NT\_EDITED.nex. For heavy-chain alignments, make sure there are no gaps within the V and within the J portion of the GERMLINE sequence (the ones that are concatenated to produce the alignment). There should be a large gap between them (corresponding to the D gene from the full VDJ sequences from the clones). We removed mostly-gap regions at the edges of the alignment.

3 - Run add\_dates\_<CLONE>.py scripts to execute python scripts that will add sequence sampling times from the original fasta files to the edited alignment files and output the final alignment files (<CLONE>\_final\_alignment.nex>)

##Deleted sequences##

During the manual edition step, we deleted the following sequences

**Alignment of VRC26 light chain**: removed sequences that contained frameshifts in the MACSE alignment	. Those are probably incomplete reads with a few "extra" or missing nucleotides at the beginning of the sequence. I removed those from the MACSE alignment and did not align them again. Removed sequences (9 out of 472 sequences, approx. 1.9%):

KJ134566
KJ134527
KJ134602
KJ134521
KJ134482
KJ134504
KJ134594
KJ134846
KJ134750

**VRC01_19** - removed sequence with frameshift indels in MACSE alignment (KP841240) Removed all-gap sites