./nextflow PlasmidSeq/pipelines/plasmid_assembly.nf \
	   --fast5 ./fast5/ \
	   --guppy_model dna_r9.4.1_450bps_hac.cfg \
	   --gpu_slot cuda:0 \
	   --tn5proj /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa.prj \
	   --bcmat /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/bc.mat \
	   -with-report report.html \
	   --tn5ref /home/f002sd4/plasmid_seq/FEX_099/demultiplexing/barcode_tn5_base.fa \
	   -resume
