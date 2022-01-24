<img align="right" src="https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true">

`CircuitSeq` is a tool to assemble and analyze plasmids from Nanopore ONT long-read sequencing. This is achieved by first basecalling the raw fast5 data and demultiplexing the sample barcodes using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), assembling with [Flye](https://github.com/fenderglass/Flye/) [Miniasm](https://github.com/lh3/miniasm), and then subsequent rounds of polishing with [Medaka](https://github.com/nanoporetech/medaka). 

## Setup CircuitSeq
Clone the repository with `git clone https://github.com/mckennalab/Circuitseq/` in the directory where you wish to perform data analysis. 

Prepare a samplesheet, an example sample sheet is provided in this repository under `example_sample_sheet.tsv`. This has to be a tab-delimited file with the following headers: `Position`, `SampleID`, `Reference`.  These correspond to the number of the barcode (e.g 1-96), the sampleID which can be a plasmid name or a alphanumeric code (avoid special characters and spaces), and optionally you can provide the location where the reference file, as a fasta, is located. This will allow the CircuitSeq pipeline to do quality assessment on the assembly. 

Modify the corresponding run file found under the pipeline directory (`runCircuitSeq.sh`). Parameter information:
```
Gotta figure out what will be removed/added, thoughts: combine the files into one run file (runCircuitSeq.sh) that can change guppy model, medaka model, barcodes

similar but different note, try to set it up so that if there is no fasta it uses a dummy pUC19 fasta reference 
```
Install [Nextflow](https://www.nextflow.io/) using the commands below:
```
#Check that you have hava 8 or later:
java -version 

#download nextflow
curl -s https://get.nextflow.io | bash 

#test that it installed correctly
./nextflow run hello 
```

Finally, once you have modified the files mentioned above, you can run the pipeline by running:
```
bash ./pipelines/runCircuitSet.sh
```

# TODO: Docker section

## Requirements
This pipeline was tested on Ubuntu 20.04. Thanks to the Docker container most dependencies are taken care of, unfortunately running Docker on some HPC nodes requires permissions which you may have to set up with your institution. 

One of the limiting factors in nanopore sequencing is the heavy GPU usage of basecalling data. It is highly recommended to run this pipeline on a computer with a powerful GPU. If your HPC node does not have GPUs it is possible to basecall the data locally and provide the basecalled reads as input to the CircuitSeq pipeline. 
```
how?
```

## Works in progress
1) Updating the pipeline to work with guppy v6+ and Q20 models, this requires some significant revamp of docker/etc which we will tackle once the Q20 stuff is available to generate test data. Apologies for the delay and thank you for your understanding.  
2) Upload a protocol file

## Our publication
For more information please refer to our publication (link). 
