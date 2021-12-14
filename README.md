![alt text](https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logosmall.png?raw=true)

`CircuitSeq` is a tool to assemble and analyze plasmids from Nanopore ONT long-read sequencing. This is achieved by first basecalling the raw fast5 data and demultiplexing the sample barcodes using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), assembling with [Miniasm](https://github.com/lh3/miniasm), and then subsequent rounds of polishing with [Racon](https://github.com/isovic/racon), [Medaka](https://github.com/nanoporetech/medaka), and [NextPolish](https://github.com/Nextomics/NextPolish). 

To run CircuitSeq:
1) Clone the repository in the directory where you wish to perform data analysis. CircuitSeq uses a docker container to simplify dependency installation, this works best if all the external files required (sample sheet, references, etc) are within the same directory where you will run CircuitSeq. 
2) Save the plasmid reference sequences in a subdirectory, these should be in the fasta format. This sequence is used for assembly quality assessment, the assembly will be made regardless of whether a reference sequence exists. If you do not have sequences for all plasmids we recommend making a dummy fasta that can be linked to each plasmid to avoid downstream quality assessment steps from interrupting the pipeline. 
3) Prepare a samplesheet, this has to be a tab-delimited file with the following headers: `Position`, `SampleID`, `Reference`. These correspond to the position on the barcode plate (i.e. the barcode ID), the sampleID which can be a plasmid name or a alphanumeric code, and the location where the reference file, as a fasta, is located. 
4) Modify the corresponding run file found under the pipeline directory (`run_941_sample_sheet.sh` for R9.4, `run_103_sample_sheet.sh` for R10.3). Items that need to be changed include: the location of the `sample sheet`, and the location of the `fast5 directory`. 
5) Install [Nextflow](https://www.nextflow.io/) using the commands below:
```
#Check that you have hava 8 or later:
java -version 

#download nextflow
curl -s https://get.nextflow.io | bash 

#test that it installed correctly
./nextflow run hello 
```

Finally, running the pipeline is as easy as:
```
bash ./pipelines/run_941_sample_sheet
```
















# OLD info: Random scripts related to Nanopore Plasmid Sequencing 

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









