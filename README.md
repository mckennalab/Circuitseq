![alt text](https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true)

`CircuitSeq` is a tool to assemble and analyze plasmids from Nanopore ONT long-read sequencing. This is achieved by first basecalling the raw fast5 data and demultiplexing the sample barcodes using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), assembling with [Miniasm](https://github.com/lh3/miniasm), and then subsequent rounds of polishing with [Racon](https://github.com/isovic/racon), [Medaka](https://github.com/nanoporetech/medaka), and [NextPolish](https://github.com/Nextomics/NextPolish). 

Setup CircuitSeq:
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
bash ./pipelines/run_941_sample_sheet.sh
```

TODO: Docker section
