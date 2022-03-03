
# Circuit-seq
<img align="right" src="https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true">

`CircuitSeq` is a pipeline to assemble and analyze plasmids from Nanopore long-read sequencing. We use a number of tools; reads are first basecalled and demultiplexed using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), and ultimately assembling with [Flye](https://github.com/fenderglass/Flye/) [Miniasm](https://github.com/lh3/miniasm), These assemblies are then polished with [Medaka](https://github.com/nanoporetech/medaka). 

# Setup 

Things you'll need:

- fast5 folder location from your Oxford Nanopore run
- Singularity (preferred) or Docker (fine, just requires admin permissions), configured with access for your username
- A server with a GPU with CUDA correctly setup (currently CUDA 11.2 is best to match TensorFlow support)

The computational pipeline starts from a directory of raw Nanopore fast5 files. Circuit-seq uses the Nextflow pipeline engine to move data through each step and output assemblies, plasmid assessments, and other information about each plasmid. 


### Install [Nextflow](https://www.nextflow.io/)

Directions are on their website: https://www.nextflow.io/, no administrator permissions needed

### Singularity setup

The tools used by Circuit-seq are often complex to install and have many dependencies. To make this easier we've packaged all the tools into a single [Docker container](https://hub.docker.com/repository/docker/aaronmck/plasmidassembly) which can be used by either the Singularity or Docker container engines. *You'll need to have either Singularity (with [user control of binds](https://singularity-admindoc.readthedocs.io/en/latest/the_singularity_config_file.html#user-bind-control-boolean-default-yes) or Docker running on the system you want to run Circuit-seq*. Unfortunately running Docker on some HPC nodes requires permissions which you may have to set up with your institution; Singularity is generally a better option on shared systems.

It's easiest to pre-download the Singularity container into the current directory (and remove it when you're done). This is done with the following command:

```
 singularity pull plasmidassembly.sif docker://aaronmck/plasmidassembly:v1_gpv5
```

This will create a file called _plasmidassembly.sif_ in the current working directory with the fully packaged Singularity container. 

### Circuit-seq setup

1. Clone the Circuit-seq repository into the directory where you'll perform data analysis:

```
git clone https://github.com/mckennalab/Circuitseq/
``` 

2. Prepare a sample sheet. An example sample sheet is provided in this repository in the [example directory](https://github.com/mckennalab/Circuitseq/tree/main/pipelines/examples). This is a tab-delimited file with the following headers: `Position`, `SampleID`, `Reference`.
  - `Position`: the number of the barcode well you used for this sample (e.g 1-96) 
  - `SampleID`: the sampleID, which can be a plasmid name or a alphanumeric code (_no_ special characters or spaces)
  - `Reference`: (optional) you can provide the location where the known fasta reference is located. Due to some Nextflow weirdness, this can't be directly in the run directory (but a subdirectory like ./references/ is fine). If you have a reference this is worth setting up, it will allow the Circuit-seq pipeline to do quality assessment on the assembly and give you aligned BAM files even when the assembly fails

3. Create a shell script and nextflow config file to run your data. You need to point to your fast5 location as well as other parameters that match you system. Here's what an example shell script looks like for running with Singularity. You'll need need to replace the bracketed values with the correct paths for your setup:

```
./nextflow  <path_to_github_checkout_of_pipelines>/pipelines/CircuitSeq.nf \
           --GPU ON \
           -c nextflow.config \
           --with-singularity \
           --samplesheet <path_to/your_sample_sheet.txt> \
           --barcodes /plasmidseq/barcodes/v2/ \
           --fast5 <full_path_to_fast5_directory> \
           --guppy_model dna_r9.4.1_450bps_sup.cfg \
           --medaka_model r941_min_sup_g507 \
           --gpu_slot cuda:0 \
           --barcode_min_score 65 \
           -resume

```

And an example nextflow.config file, filling in _your_singularity_sif_file_path_here_ with your path to the _.sif_ from the ```singularity pull``` command :

```
params.quality_control_processes = true
params.nextpolish_cfg = "<path_to_github_checkout_of_pipelines>/pipelines/run.cfg"

singularity{
        enabled = true
        autoMounts = true
}
process {
        container = '<your_singularity_sif_file_path_here>'
  withLabel: with_gpus {
         maxForks = 1
         containerOptions = { workflow.containerEngine == "singularity" ? '--nv':
         ( workflow.containerEngine == "docker" ? '--gpus all': null )}
  }
}
```


4. Finally, once you have modified the files mentioned above, you can run the pipeline by running:
```
bash <shell_script_you_saved_just_above_here.sh>
```

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
