
# Circuit-seq
<img align="right" src="https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true">

`CircuitSeq` is a pipeline to assemble and analyze plasmids from Nanopore long-read sequencing. We use a number of tools; reads are first basecalled and demultiplexed using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), and ultimately assembling with [Flye](https://github.com/fenderglass/Flye/) [Miniasm](https://github.com/lh3/miniasm), These assemblies are then polished with [Medaka](https://github.com/nanoporetech/medaka). 

## Setup 

Things you'll need:

- fast5 folder location from your Oxford Nanopore run
- Docker, configured with access for your username
- Optionally a GPU, which speeds up basecalling and other steps

The computational pipeline starts from a directory of raw Nanopore fast5 files. Circuit-seq uses the Nextflow pipeline engine to shuffle data through each step and output assemblies, plasmid assessments, and other information about each plasmid. First you'll need to have Nextflow installed:

### Install [Nextflow](https://www.nextflow.io/) using the commands below:

- Check that you have have version 8 or later:
```java -version ```

- download nextflow:
```curl -s https://get.nextflow.io | bash ```

- test that it installed correctly
```./nextflow run hello ```

If you have issues, [Nextflow](https://www.nextflow.io/) has great community resources linked from their website.

### Singularity setup

The tools used by Circuit-seq are often complex to install and have many dependencies. To make this easier we've packaged all the tools into a single [Docker container](https://hub.docker.com/repository/docker/aaronmck/plasmidassembly) which can be used by either the Singularity or Docker container engines. *You'll need to have either Singularity (with [user control of binds](https://singularity-admindoc.readthedocs.io/en/latest/the_singularity_config_file.html#user-bind-control-boolean-default-yes) or Docker running on the system you want to run Circuit-seq*. Unfortunately running Docker on some HPC nodes requires permissions which you may have to set up with your institution; Singularity is generally a better option on shared systems.

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

## Learn more / cite our publication

For more information please refer to [our publication](https://www.biorxiv.org/content/10.1101/2022.01.25.477550v1). If you find this work useful, also consider citing our work: 

```
@article {Emiliani2022.01.25.477550,
	author = {Emiliani, Francesco E and Hsu, Ian and McKenna, Aaron},
	title = {Circuit-seq: Circular reconstruction of cut in vitro transposed plasmids using Nanopore sequencing},
	elocation-id = {2022.01.25.477550},
	year = {2022},
	doi = {10.1101/2022.01.25.477550},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Recombinant DNA is a fundamental tool in biotechnology and medicine. Validation of the resulting plasmid sequence is a critical and time-consuming step, which has been dominated for the last 35 years by Sanger sequencing. As plasmid sequences grow more complex with new DNA synthesis and cloning techniques, we need new approaches that address the corresponding validation challenges at scale. Here we prototype a high-throughput plasmid sequencing approach using DNA transposition and Oxford Nanopore sequencing. Our method, Circuit-seq, creates robust, full-length, and accurate plasmid assemblies without prior knowledge of the underlying sequence for approximately $1.50 per plasmid. We demonstrate the power of Circuit-seq across a wide range of plasmid sizes and complexities, generating accurate and contiguous plasmid maps. We then leverage our long read-data to characterize epigenetic marks and estimate plasmid contamination levels. Circuit-seq scales to large numbers of samples at a lower cost than commercial Sanger sequencing, accelerating a key step in synthetic biology, with low startup costs make it practical for individual laboratories.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2022/01/26/2022.01.25.477550},
	eprint = {https://www.biorxiv.org/content/early/2022/01/26/2022.01.25.477550.full.pdf},
	journal = {bioRxiv}
}
```
