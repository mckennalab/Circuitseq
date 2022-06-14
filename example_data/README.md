# Test runs
This is an example dataset, files, and shell script to run the pipeline directly from raw fast5s or from pre-basecalled fastqs. 

The purpose of this set is to give you a chance to test that CircuitSeq will work on your system. On our system the whole process from basecalling takes ~2 minutes. 

To run the two test cases you need to:
1. Follow the steps in our main readme to set up nextflow, singularity, and clone CircuitSeq 
        * Head to `run_with_basecalling` if you want to try running with basecalling, modify the `run_nf.sh` with the appropriate paths and run it with `bash run_nf.sh`
        * Head to `run_starting_with_fastqs` if you want to try running with pre-basecalled data, modify the `run_nf.sh` with the appropriate paths, note that here you will have to provide two new parameters, and run it with `bash run_nf.sh`

## What to expect?
We kept this simple because fast5s are large and we wanted this to run quickly. There are downsampled reads just for one plasmid ([LentiCas9-Blastidicin](https://www.addgene.org/52962/)).  
