# PlasmidSeq
Random scripts related to Nanopore Plasmid Sequencing 

Link to google drive for figures: https://drive.google.com/drive/folders/1cBPizfFrt0MZ8FUU1ypLNDUxe284hu79?usp=sharing

link to paper https://docs.google.com/document/d/16PynW7_z6Mo8nsEQhlPXQrMAuR2YQQ6ERzj-_hfyYp8/edit?ts=602ea4b9



For the second 96-well run (that worked a little better)

location of raw data:
/dartfs/rc/lab/M/McKennaLab/projects/Oxford_computer_backup/data/FEX_092/

Location of output
/dartfs/rc/lab/M/McKennaLab/projects/Plasmid_Sequencing/plasmidseq/FEX_092/

Steps involved:
1) basecall with guppy4.5.2 model = dna_r9.4.1_450bps_hac.cfg 
2) go into the demultiplexing dir this process was adapted from https://www.protocols.io/view/demultiplexing-nanopore-reads-with-last-7vmhn46?step=7
   In it you will find a barcode_tn5_base.fa with the barcode sequences (no additional context, just the barcode) which was indexed with:
   lastdb -uNEAR barcode_tn5_base.fa barcode_tn5_base.fa
   and a bc.mat substitution matrix 
   The script demultiplex_samples.sh has the following steps for demultiplexing:
   1) extracts the front and back end of each read and concatenates them using extract_ends.py
   2) using lastal you assign the barcodes
   3) then you bin them using the fastx-fetch.pl script 
   4) and finally the Rscript makes those cool graphs using the lengths_* files but maybe we need to find a way to use only the lengths file submitted by the 
   same user and limit it  to 6? plasmids per set because it gets crowded up in there 
3) the demultiplexing/demultiplexed dir now holds the demultiplexed fastq and the lengths files 
4) in the demultiplexing/porechop folder you run the porechop_trim_and_unchimera.sh script to find the nanopore adapters and trim them + extra, also if it finds them in the middle of the read it tosses the read out since that is the product of two plasmids that got ligated (we could try to then subsequently bin them correctly but really it's quite rare. 
5) in the demultiplexing/filtered dir filtering.sh runs the filtering using NanoFilt filters by "quality" currently set at 10 (uses the sequencing_summary.txt file output by the basecalling) 
6) Now we are ready to assemble!
7) in the FEX_092/assemblies dir you will find the assembly_for_iterator.sh file... this one is a little wild
   1) canu correct the reads, this reduces the number of reads this 100% speeds up the process, in the past i saw it was giving me better assemblies... but maybe we can try removing it now? 
   2) map the reads to each other to create overlaps with minimap2
   3) ***use miniasm to assemble the options -o 3000 -I 0.9 -F 0.9 were added by me because i noticed it improved assemblies. these would be good ones to be able to modify to see if we can 100% optimize the assembly process***
   4) convert the .gfa to fasta and racon polish 3 times (this is a pretty common number of racon iterations, but we could set up a modifiable "pick how many times you want to racon option" 
   5) do a final polish with medaka **here the model used is very important** r941_min_high_g360 is the correct model for the 9.4 basecalled with the latest medaka, but this could 1) change and 2) will be different for the 10.3 so this should be an option we can set in the pipeline. 
   6) pilon polish... this was added recently, i noticed that sometimes medaka caused some smallish deletions (like 5bp or so) but pilon cleans that up 
   7) pilon outputs those stupid fasta files where each line is like 50 bases and then it creates a new line, this just makes it one block of text you can copy paste easily, because i am lazy. and renames it since at this point all file names are the same for all plasmids 
   8) done! 









