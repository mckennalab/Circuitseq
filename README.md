
# Circuit-seq
<img align="right" src="https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true">

`CircuitSeq` is a pipeline to assemble and analyze plasmids from Nanopore long-read sequencing. We use a number of tools; reads are first basecalled and demultiplexed using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), and ultimately assembling with [Flye](https://github.com/fenderglass/Flye/) [Miniasm](https://github.com/lh3/miniasm), These assemblies are then polished with [Medaka](https://github.com/nanoporetech/medaka). 

## Setup 

Things you'll need:

- fast5 folder location from your Oxford Nanopore run
- Docker setup with access for your username
- Optionally a GPU, which speeds up basecalling and other steps

The computational pipeline starts from a directory of raw Nanopore fast5 files. Circuit-seq uses the Nextflow pipeline engine to shuffle data through each step and output assemblies, plasmid assessments, and other information about each plasmid. First you'll need to have Nextflow installed:

### Install [Nextflow](https://www.nextflow.io/) using the commands below:

- Check that you have hava 8 or later:
```java -version ```

- download nextflow:
```curl -s https://get.nextflow.io | bash ```

- test that it installed correctly
```./nextflow run hello ```

If you have issues, [Nextflow](https://www.nextflow.io/) has great community resources linked from their website.

### Docker setup

The tools used by Circuit-seq are often complex to install and have many dependencies. To make this easier we've packaged all the tools into a single [Docker container](https://hub.docker.com/repository/docker/aaronmck/plasmidassembly). *You'll need to have Docker running on the system you want to run Circuit-seq*. Unfortunately running Docker on some HPC nodes requires permissions which you may have to set up with your institution. We're also working on testing support for [Singularity](https://sylabs.io/guides/3.5/user-guide/introduction.html) which is generally easier to run on distributed systems. 

### Setup Circuit-seq

- Clone the repository with ```git clone https://github.com/mckennalab/Circuitseq/``` in the directory where you wish to perform data analysis. 

- Prepare a samplesheet. An example sample sheet is provided in this repository under `example_sample_sheet.tsv`. This is a tab-delimited file with the following headers: `Position`, `SampleID`, `Reference`.  
  - `Position`: the number of the barcode (e.g 1-96) 
  - `SampleID`: the sampleID which can be a plasmid name or a alphanumeric code (avoid special characters and spaces)
  - `Reference`: (optional) you can provide the location where the known fasta reference is located. This will allow the Circuit-seq pipeline to do quality assessment on the assembly

- Copy and modify the example run file found under the pipeline directory (`runCircuitSeq.sh`) to point to your fast5 location as well as other parameters that match you system.

- Finally, once you have modified the files mentioned above, you can run the pipeline by running:
```
bash ./pipelines/runCircuitSet.sh
```

## Works in progress / TODOs:
- See the Github issues tab above for future computational improvements
- Support for guppy v6+ and Q20 models; this will greatly improve assemblies
- Cleaner packaging for AWS or other cloud computing platforms
- Visualization of the final assemblies

## Our publication
For more information please refer to our publication (link). 
