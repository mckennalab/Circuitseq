# replace <singularity.sif> with the absolute path to your singularity container
# also be sure to run this script inside the example_data/run_with_basecalling directory
NXF_VER=21.10.6 ./nextflow run ../../pipelines/CircuitSeq.nf \
           --GPU ON \
           -c ../../pipelines/nextflow.config \
           -with-singularity <singularity.sif> \
           --samplesheet runnable_example_samplesheet.tsv \
           --use_existing_basecalls false \
           --fast5 ../fast5/ \
           --basecalling_dir "" \
           --base_calling_summary_file "" \
           --barcodes /plasmidseq/barcodes/v2/ \
           --barcode_kit "MY-CUSTOM-BARCODES" \
           --guppy_model dna_r9.4.1_450bps_sup.cfg \
           --medaka_model r1041_e82_400bps_hac_g615 \
           --gpu_slot cuda:0 \
           --barcode_min_score 65 \
           --quality_control_processes true \
           -resume
