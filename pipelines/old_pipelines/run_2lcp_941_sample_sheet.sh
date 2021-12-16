
./nextflow ./pipelines/plasmid_assembly_2fixed_duplications.nf \
	   --samplesheet ./fex_120_sample_sheet.tsv \
	   --barcodes /analysis/2021_03_12_nanopore_stuff/guppy_plasmid_barcodes_27mer \
	   --fast5 /dartfs2/rc/lab/M/McKennaLab/projects/Oxford_computer_backup/data/FEX_120/96_v2/20210803_1913_MN34875_AHC239_bdcd3859/fast5/ \
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
