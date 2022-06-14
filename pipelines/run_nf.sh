#It is safest to use absolute paths  
NXF_VER=21.10.6 nextflow run ./Experimental_Circuitseq/pipelines/CircuitSeq.nf \
           --GPU ON \
           -c nextflow.config \
           -with-singularity <path to .sif file> \
           --samplesheet <path to sample_sheet.tsv> \
           --fast5 <path to fast5 directory> \
           --barcodes /plasmidseq/barcodes/v2/ \
           --guppy_model dna_r9.4.1_450bps_sup.cfg \
           --medaka_model r941_min_sup_g507 \
           --gpu_slot cuda:0 \
           --barcode_min_score 65 \
           --quality_control_processes true \
           --use_existing_basecalls false \
           -resume

#if you are running from previously basecalled data and want to skip basecalling change "use_existing_basecalls" to true 
#add the following parameters to this file:
#path to fastq directory
#--basecalling_dir <path_to_fastq_dir>
#path to sequencing_summary file
#--base_calling_summary_file <path_to_summary.txt>
