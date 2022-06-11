# Test runs
This is an example dataset, run scripts, and configuration files to run the pipeline directly from raw fast5s or from pre-basecalled fastqs. 

The purpose of this set is to give you a chance to test that CircuitSeq will work on your system, but also for me to explain a quirk of `Nextflow` that has caused people problems. Namely that the `nextflow.config` file: 
1. Requires **absolute paths**, and wont accept relative paths. 
2. Does not accept empty parameters. So if you are running the pipeline with basecalling you need to **comment out the two parameters used for running with fastqs**. 
3. Paths in the `run_nf.sh` file can be relative. 

To run the two test cases you need to:
1. Follow the steps in our main readme to set up CircuitSeq
2. Move into either the `run_with_basecalling` or `run_starting_with_fastqs` 
  * If you want to try running with basecalling, you can run the pipeline directly with `bash run_nf.sh`.  The `nextflow.config` file is already set up correctly with `params.use_existing_basecalls` set to `false` and `params.basecalling_dir` and `params.base_calling_summary_file` commented out. Because the paths in the `run_nf.sh` are relative those also dont need to be changed.
  * If you are running with basecalled data (i.e. fastqs) you will need to change the paths in the `run_starting_with_fastqs/nextflow.config` for the fastqs and the sequencing_summary.txt to include where you cloned this repository. 

## What to expect?
We kept this simple because fast5s are large and we wanted this to run quickly. There are downsampled reads just for one plasmid ([LentiCas9-Blastidicin](https://www.addgene.org/52962/)).  
