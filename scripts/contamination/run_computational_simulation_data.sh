
../.././nextflow ../../pipelines/plasmid_assembly_from_filtered_reads.nf \
	   --GPU ON \
	   -c nextflow.config \
	   --with-docker \
	   --basecalling_summary_file sequencing_summary.txt \
	   --samplesheet computational_96_well_mixture_addgene.txt \
	   --barcodes /plasmidseq/barcodes/v1/ \
	   --fast5 /analysis/2021_07_12_PlasmidSeq_publication_96xr941/fast5/ \
	   --guppy_model dna_r9.4.1_450bps_sup.cfg \
	   --medaka_model r941_min_sup_g507 \
	   --gpu_slot cuda:0 \
	   --tn5proj /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa.prj \
	   --bcmat /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/bc.mat \
	   -with-report report.html \
	   --tn5ref /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa \
	   --barcode_min_score 65 \
	   --nanopolish /analysis/2021_05_06_nanopore_pipeline/NextPolish/nextPolish \
	   --nanopolish_run_config /analysis/2021_05_06_nanopore_pipeline/bin/run.cfg \
	   -resume
# -with-docker 065c7ae93e83 \
# dna_r9.4.1_450bps_sup.cfg \
#
#           --bcmat /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/bc.mat \
# /analysis/2021_08_26_PlasmidSeq_paper/fast5
