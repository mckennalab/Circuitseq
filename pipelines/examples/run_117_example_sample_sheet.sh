
./nextflow <path_to_github_checkout_of_pipelines>/pipelines/CircuitSeqFlyMed.nf \
	   --GPU ON \
	   -c nextflow.config \
	   --with-docker \
	   --samplesheet ./example_sample_sheet.tsv \
	   --barcodes /plasmidseq/barcodes/v2/ \
	   --fast5 <path_to_fast5_file_directory> \
	   --guppy_model dna_r9.4.1_450bps_sup.cfg \
	   --medaka_model r941_min_sup_g507 \
	   --gpu_slot cuda:0 \
	   -with-report report.html \
	   --barcode_min_score 65 \
	   -resume
