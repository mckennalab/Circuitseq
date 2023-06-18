
# Circuit-seq
<img align="right" src="https://github.com/mckennalab/Circuitseq/blob/main/circuitSeq_logo_red.png?raw=true">

`CircuitSeq` is a pipeline to assemble and analyze plasmids from Nanopore long-read sequencing. We use a number of tools; reads are first basecalled and demultiplexed using [Guppy](https://nanoporetech.com/), filtering out chimeric reads and short reads (with [Porechop](https://github.com/rrwick/Porechop) and [NanoFilt](https://github.com/wdecoster/nanofilt)), correcting with the reads with [Canu](https://github.com/marbl/canu), and ultimately assembling with [Flye](https://github.com/fenderglass/Flye/) [Miniasm](https://github.com/lh3/miniasm), These assemblies are then polished with [Medaka](https://github.com/nanoporetech/medaka). 

# Setup 

Things you'll need if running from fast5s:

- location of your fast5 from your Oxford Nanopore run
- A server with a GPU with CUDA correctly setup (currently CUDA 11.2 is best to match TensorFlow support)
- Singularity (preferred) or Docker (fine, just requires admin permissions), configured with access for your username


Things you'll need if running with basecalled data:
- location of your fastq directory 
- location of the sequencing_summary.txt file (usually found in the basecalling dir)
- Singularity (preferred) or Docker (fine, just requires admin permissions), configured with access for your username


The computational pipeline starts from a directory of raw Nanopore fast5 files, you can also skip the basecalling step if you already have basecalled your data. Circuit-seq uses the Nextflow pipeline engine to move data through each step and output assemblies, plasmid assessments, and other information about each plasmid. 


### Install [Nextflow](https://www.nextflow.io/)

Directions are on their website: https://www.nextflow.io/, no administrator permissions needed
Ideally put nextflow in your PATH

### Singularity setup

The tools used by Circuit-seq are often complex to install and have many dependencies. To make this easier we've packaged all the tools into a single [Docker container](https://hub.docker.com/repository/docker/aaronmck/plasmidassembly) which can be used by either the Singularity or Docker container engines. *You'll need to have either Singularity (with [user control of binds](https://singularity-admindoc.readthedocs.io/en/latest/the_singularity_config_file.html#user-bind-control-boolean-default-yes) or Docker running on the system you want to run Circuit-seq*. Unfortunately running Docker on some HPC nodes requires permissions which you may have to set up with your institution; Singularity is generally a better option on shared systems.

It's easiest to pre-download the Singularity container into the current directory (and remove it when you're done). This is done with the following command:

```
 singularity pull plasmidassembly.sif docker://aaronmck/plasmidassembly:1_0_1
```

This will create a file called _plasmidassembly.sif_ in the current working directory with the fully packaged Singularity container. 

Older versions of singularity (v2 and below) can be downloaded with:
```
 singularity build plasmidassembly.sif docker://aaronmck/plasmidassembly:1_0_1
```

### Circuit-seq setup

1. Clone the Circuit-seq repository into the directory where you'll perform data analysis:

```
git clone https://github.com/mckennalab/Circuitseq/
``` 

2. Prepare a sample sheet. An example sample sheet is provided in this repository in the [example directory](https://github.com/mckennalab/Circuitseq/tree/main/example_data/example_samplesheet.tsv). This is a tab-delimited file with the following headers: `position`, `sample`, `reference`.
  	- `position`: the number of the barcode well you used for this sample (e.g 01-96) 
  	- `sample`: the sampleID, which can be a plasmid name or a alphanumeric code (_no_ special characters or spaces)
  	- `reference`: you can provide the location where the known fasta reference is located. Due to some Nextflow weirdness, this can't be directly in the run directory (but a subdirectory like ./references/ is fine). If you have a reference, this is worth setting up, it will allow the Circuit-seq pipeline to do quality assessment on the assembly and give you aligned BAM files even when the assembly fails. If you don't have a reference simple fill in this column with `NA`.

3. Create a copy of the `run_nf.sh` shell script found in either one of the [example data](https://github.com/mckennalab/Circuitseq/tree/main/example_data) directory and modify the following parameters:
    - Path to nextflow if it is not in your path
    - Path to the CircuitSeq.nf pipeline file found in /pipelines 
    - Path to the nextflow.config file found in /pipelines
    - Path to your samplesheet 
    - Choose if you are running from fast5 or fastq by changing use_existing_basecalls to false or true, respectively 
    - If running from fast5 provide a path to `--fast5` directory and leave `--basecalling_dir` and `--base_calling_summary_file` as `""`
    - If running from fastq leave `--fast5`as `""` and provide fastq directory for`--basecalling_dir` and sequencing_summary.txt file for`--base_calling_summary_file`


#To use nanopore barcoding kits you can change: 
#--barcodes to /plasmidseq/barcodes/nanopore_official/
#--barcode_kit to the name of the barcode kit you used (this is with guppy v5.0.16 names which can be found in our barcodes/nanopore_official directory on the github. 
#our --barcode_min_score was set based on our demultiplexing data to achieve best sensititivy/specificity, if you find that using nanopore barcodes you are losing too many reads or getting too much noise you can change this parameter accordingly. 

```

4. Finally, once you have modified the files mentioned above, you can run the pipeline by running:
```
bash <shell_script_you_modified_just_above_here.sh>
```

### Barcodes and custom barcodes
The nanopore official barcodes are already included in the docker container. If you use the nanopore kits you need to set the following parameters:
```
 --barcodes /plasmidseq/barcodes/nanopore_official/ \
 --barcode_kit "kit-name-in-quotes" \
```
The kit name is the one you would use for normal guppy demulitplexing, eg: `"SQK-RBK110-96"`

If you would like to use your own barcodes you can by mounting your custom files in the singularity container by adding `--bind /source/directory:/location/in/container` to the run script and then changing `--barcodes` to point to where you mounted them.  e.g.:

```
           --barcodes /mnt/barcodes/ \
           --barcode_kit "MY-NEW-CUSTOM-BARCODES" \
```

### Running 10.4.1 Flongles

Circuitseq was developed with the 9.4 chemistry but supports the current 10.4.1 Flongles. You'll just need to change the base calling models appropriately:

```
--guppy_model dna_r9.4.1_450bps_sup.cfg \
--medaka_model r1041_e82_400bps_hac_g615 \
```

### Test data
If you want to test out the pipeline we have added a downsampled fast5 and fastqs in the `example_data` directory. There you will find an explanation of how to run each example. 

## pLannotate your plasmids
By popular demand we are trying to implement a quick way of annotating your plasmids using Matt McGuffie's [pLannotate](https://github.com/barricklab/pLannotate). If you are interested please go to the [plannotate](https://github.com/mckennalab/Circuitseq/edit/main/plannotate) directory and follow the readme. 


## Learn more / cite our publication

For more information please refer to [our publication](https://pubs.acs.org/doi/10.1021/acssynbio.2c00126). If you find this work useful, also consider citing our work: 

```
@article{doi:10.1021/acssynbio.2c00126,
author = {Emiliani, Francesco E. and Hsu, Ian and McKenna, Aaron},
title = {Multiplexed Assembly and Annotation of Synthetic Biology Constructs Using Long-Read Nanopore Sequencing},
journal = {ACS Synthetic Biology},
volume = {0},
number = {0},
pages = {null},
year = {0},
doi = {10.1021/acssynbio.2c00126},
    note ={PMID: 35695379},

URL = {
        https://doi.org/10.1021/acssynbio.2c00126

}
```
