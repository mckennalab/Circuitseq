
./nextflow ./pipelines/plasmid_assembly_fixed_duplications.nf \
	   --samplesheet ./fex_119_sample_sheet.tsv \
	   --barcodes /analysis/2021_03_12_nanopore_stuff/guppy_plasmid_barcodes_27mer \
	   --fast5 /home/f002sd4/plasmid_seq/FEX_119/FEX_119/New_barcodes/20210801_2032_MN34875_AHC194_91e39b7f/fast5 \
	   --guppy_model dna_r9.4.1_450bps_sup.cfg \
	   --medaka_model r941_min_sup_g507 \
	   --gpu_slot cuda:0 \
	   --tn5proj /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa.prj \
	   --bcmat /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/bc.mat \
	   -with-report report.html \
	   --tn5ref /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa \
	   --barcode_min_score 65 \
	   --nanopolish /analysis/2021_05_06_nanopore_pipeline/NextPolish/nextPolish \
           --nanopolish_run_config /analysis/2021_03_12_nanopore_stuff/run.cfg \
	   -resume
