#!/usr/bin/env nextflow

/*
========================================================================================
                         mckenna_lab/plasmid_seq
========================================================================================
 Process plasmids sequenced on Oxford Nanopore into assembled maps of each sequence
----------------------------------------------------------------------------------------
*/


input_fastq5_path  = Channel.fromPath(params.fast5)
input_tn5ref      = Channel.fromPath(params.tn5ref)
input_tn5proj_base = Channel.fromPath("${params.tn5ref}.*")

// input_tn5proj_base.subscribe { println "value: $it" }
input_bc_mat     = Channel.fromPath(params.bcmat)
results_path = "results"
minimum_file_size = 5000
minimum_fastq_size = 50000


/*
 * Read in the sample table and convert entries into a channel of sample information 
 */

Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.position.padLeft(2,'0'), it.sample, file(it.reference) ) }
        .set { sample_table }

//sample_table_print.subscribe { println "\n[samples_R1_R2] ${it}\n"}


/*
 * Basecalling using Guppy
 *
 * The main output is a channel input_fastq5_path which contains a merged fastq.gz file of all of passing reads
 */
process GuppyBaseCalling {
    publishDir "$results_path/guppy"

    input:
    path input_path from input_fastq5_path

    output:
    path "basecalling/pass/" into basecalled_directory
    path "basecalling/sequencing_summary.txt" into basecalling_summary   // the fastq output file path

    script:
        
    """
    guppy_basecaller \\
        --input_path $input_path \\
        --save_path ./basecalling \\
        -c $params.guppy_model \\
        -x $params.gpu_slot \\
        --records_per_fastq 0 \\
        --compress_fastq
        
    guppy_basecaller --version &> v_guppy.txt
    """
}

process GuppyDemultiplex {
    publishDir "$results_path/guppy_demultiplex"

    input:
    path basecalled from basecalled_directory

    output:
    path "saved_data/barcoding_summary.txt" into barcoding_split_summary   // the fastq output file path
    path "saved_data/barcode*/*.fastq.gz" into fastq_gz_split_files mode flatten

    script:
        
    """
    guppy_barcoder \\
        --input_path ${basecalled} \\
        --save_path saved_data \\
        --data_path ${params.barcodes} \\
        --barcode_kits MY-CUSTOM-BARCODES \\
        --front_window_size 120 \\
        --min_score_mask 30 \\
        -q 1000000 \\
        --min_score 65 \\
        --compress_fastq
    """
}

process LengthFilter {
    publishDir "$results_path/length_filter"
    
    input:
    tuple datasetID, file(fastq) from fastq_gz_split_files.filter(){ it.countFastq() > 5} .map { file -> tuple( (file.toString().split("barcode"))[1][0..1], file) }

    output:
    tuple datasetID,"${datasetID}_length_filtered.fq.gz" into length_output
    script:
        
    """
    filter_by_length.py \
    ${fastq} \
    ${datasetID}_length_filtered.fq.gz
    """
}



process Porechop {
    publishDir "$results_path/porechop"
    
    input:
    tuple val(datasetID), file(fastq) from length_output
    output:
    tuple datasetID,"${datasetID}_porechop.fq.gz" into porechop_output
    script:
        
    """
    /home/f002sd4/plasmid_seq/Porechop_for_plasmidseq/porechop-runner.py -i ${fastq} -o ${datasetID}_porechop.fq.gz \
                                                            --format fastq.gz \
                                                            --end_threshold  50 --extra_end_trim 10 \
                                                            --discard_middle --middle_threshold 80
    """
}


// this is very nextflow -- we're going to consume the summary file many times (96), and we 
// need to convert it into a value channel which can be reused again and again, instead of a 
// queue channel. 
basecalling_summary_file = basecalling_summary.first()

// TODO: nanofilt needs to be installed
process FilterReads {
    publishDir "$results_path/filter_reads"
    
    input:
    tuple val(datasetID), file(tofilter) from porechop_output
    path summary from basecalling_summary_file
    
    output:
    tuple val(datasetID), file("${datasetID}_filtered.fq.gz") into filtered_reads, filtered_reads_racon, filtered_reads_racon2, filtered_reads_racon3, filtered_reads_medaka, filtered_reads_pilon, filtered_reads_minimap2, filtered_reads_minimap3

    script:
        
    """
    gunzip -c ${tofilter} | NanoFilt --quality 10 --length 2500 --summary ${summary} | gzip > ${datasetID}_filtered.fq.gz
    """
}

// TODO: fix path to canu
// TODO: Canu really struggles with low read counts, right now we're filtering out anything with a file size less than 100Kb
process CanuCorrect {
    publishDir "$results_path/canu"
    
    input:
    tuple val(datasetID), file(to_correct) from filtered_reads //.filter{ it.get(1).size()>minimum_file_size}
    
    output:
    tuple val(datasetID), file("${datasetID}_canu_correct/reads.correctedReads.fasta.gz") into canu_corrected_minimap, canu_corrected_convert
    
    script:

    """
    mkdir ${datasetID}_canu_correct
    /analysis/2021_05_06_nanopore_pipeline/canu-2.1.1/bin/canu -correct \
     -p reads -d ${datasetID}_canu_correct \
     genomeSize=8k \
     stopOnLowCoverage=2 minInputCoverage=2 \
     -nanopore ${to_correct} 
    
    
    """
}

// TODO: install miniasm in a better way
process Miniasm {
    publishDir "$results_path/miniasm"
    
    input:
    tuple val(datasetID), path(corrected_reads) from canu_corrected_minimap

    output:
    tuple val(datasetID), path("${datasetID}_overlaps.gfa") into miniasm_overlaps
    
    script:
    """
    minimap2 -x ava-ont \
         ${corrected_reads} \
         ${corrected_reads} > ${datasetID}_overlaps.paf

    miniasm \
    -o 3000 -I 0.9 -F 0.9 \
    -f ${corrected_reads} ${datasetID}_overlaps.paf > ${datasetID}_overlaps.gfa
    """
}

read_phased = canu_corrected_convert.phase(miniasm_overlaps)

process ConvertGraph {
    errorStrategy 'finish'
    publishDir "$results_path/convert_graph"
    
    input:
    val phased from read_phased.filter(){it.get(1).get(1).countLines() > 1} //[val(datasetID),path(corrected_reads)],[val(datasetID),val(miniasm_overlap)]] .filter{ tag, file -> !file.isEmpty() }
    
    output:
    tuple val(str_name), path("${str_name}_contigs.fasta") into fasta_graph

    script:
    str_name = phased.get(0).get(0)
    $/
    #!/usr/bin/env python
    import sys
 
    ov = open("${phased.get(1).get(1)}")
    ln = ov.readline()
    spl = ln.split()
    n_prop = sum([1 if x == 'N' else 0 for x in spl[2]])/len(spl[2])

    if n_prop < 0.5 and len(spl[2]) > 100:
        output = open("${phased.get(0).get(0)}_contigs.fasta","w")
        output.write(">" + spl[1] + '\n')
        output.write(spl[2] + '\n')
        output.close()
    else:
        output = open("${phased.get(0).get(0)}_contigs.fasta","w")
        output.close()
    /$
}

read_phased_racon = filtered_reads_racon.phase(fasta_graph)


// TODO: install racon the right way
process RaconPolish {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
   
    input:
    val tuple_pack from read_phased_racon.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon1.fasta") into racon_corrected
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -m 8 -x -6 -g -8 -w 500 -t 8 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon1.fasta
    """
}



read_phased_racon2 = filtered_reads_racon2.phase(racon_corrected)

// TODO: install racon the right way
process RaconPolish2 {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
   
    input:
    val tuple_pack from read_phased_racon2.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon2.fasta") into racon_corrected2
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -m 8 -x -6 -g -8 -w 500 -t 8 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon2.fasta
    """
}



read_phased_racon3 = filtered_reads_racon3.phase(racon_corrected2)

// TODO: install racon the right way
process RaconPolish3 {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
   
    input:
    val tuple_pack from read_phased_racon3.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon3.fasta") into racon_corrected3
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -m 8 -x -6 -g -8 -w 500 -t 8 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon3.fasta
    """
}

read_phased_medaka = filtered_reads_medaka.phase(racon_corrected3)

// TODO: install medaka_consensus the right way
// TODO: parameter for the model
/* TODO: requires:
bcftools   Not found  1.11       False
bgzip      Not found  1.11       False
minimap2   2.17       2.11       True
samtools   1.10       1.11       False
tabix      Not found  1.11       False
*/
process MedakaConsensus {
    errorStrategy 'finish'
    publishDir "$results_path/medaka_consensus"
   
    input:
    val tuple_pack from read_phased_medaka.filter{ it.get(1).get(1).countFasta()>=1} 
    
    output:
    set str_name, file("${str_name}_racon_medaka/consensus.fasta") into racon_medaka_consensus
    
    script:
    str_name = tuple_pack.get(0).get(0)
    """
    medaka_consensus -i ${tuple_pack.get(0).get(1)} -d ${tuple_pack.get(1).get(1)} -o ${str_name}_racon_medaka -m r941_min_high_g360

    """
}

process Pilon {
    errorStrategy 'finish'
    publishDir "$results_path/pilon"
   
    input:
    set datasetID, file(filtered_reads) from filtered_reads_pilon
    set dataset_graph, file(medaka_consensus) from racon_medaka_consensus

    output:
    set datasetID, file("${dataset_graph}_consensus.fasta") into pilon_consensus
    
    script:

    """
    minimap2 -ax map-ont ${medaka_consensus} ${filtered_reads} | samtools view - -Sb | samtools sort - -@8 -o mapping.sorted.bam
    samtools index mapping.sorted.bam

    java -Xmx32G -jar /analysis/2021_05_06_nanopore_pipeline/bin/pilon-1.24.jar \
        --genome ${medaka_consensus} --fix all --changes \
        --unpaired mapping.sorted.bam --output ${dataset_graph}_consensus

    """
}
//pilon_consensus.view { "value: $it" }
//sample_table.view { "sample_table: $it" }
read_phased_pilon = pilon_consensus.phase(sample_table)

process MarsCorrection {
    errorStrategy 'finish'
    publishDir "$results_path/mars"
   
    input:
    val tuple_pack from read_phased_pilon
    //set samplePosition, sample_ID, reads1, reads2  from sample_table

    
    output:
    set str_name, path("${str_name}_rotated.fasta") into mars_rotated
    
    script:
    str_name = tuple_pack.get(1).get(1)

    """

    cat ${tuple_pack.get(1).get(2)} ${tuple_pack.get(0).get(1)} > merged_reference.fa
    
    /analysis/2021_04_22_circular_alignment/MARS/mars -a DNA -i merged_reference.fa -o ${tuple_pack.get(1).get(1)}_rotated.fasta -m 1

    """
}