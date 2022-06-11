NXF_VER=21.10.6 nextflow run /analysis/Francesco/Plasmid_Sequencing/2022_06_07_CircuitSeq_test_data/Circuitseq/pipelines/CircuitSeq.nf \
           --GPU ON \
           -c nextflow.config \
           --with-singularity \
           --samplesheet /analysis/Francesco/Plasmid_Sequencing/2022_06_07_CircuitSeq_test_data/Circuitseq/example_data/example_tearsheet.tsv \
           --barcodes /plasmidseq/barcodes/v2/ \
           --fast5 /analysis/Francesco/Plasmid_Sequencing/2022_06_07_CircuitSeq_test_data/Circuitseq/example_data/fast5 \
           --guppy_model dna_r9.4.1_450bps_sup.cfg \
           --medaka_model r941_min_sup_g507 \
           --gpu_slot cuda:0 \
           --barcode_min_score 65 \
           -resume
