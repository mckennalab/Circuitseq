#!/usr/bin/env nextflow

/*
conda create -n medaka -c conda-forge -c bioconda medaka
conda activate medaka
conda install -c bioconda nanofilt
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
        .map{ tuple(it.position.padLeft(2,'0'), it.sample, file(it.reference), it.sangers ) }
        .into{sample_table; sample_table_assessment; sample_table_print}

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
    path "saved_data/barcode**/**.fastq.gz" into fastq_gz_split_files

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


fastq_gz_split_files.flatten().filter(){ it.countFastq() > 1 && !(it.getParent().contains("unclassified"))}.into{ guppy_demulti; guppy_demulti_print }
//guppy_demulti_print.flatten().view()

process LengthFilter {
    publishDir "$results_path/length_filter"
    
    input:
    tuple datasetID, path(fastq) from guppy_demulti.map { file -> tuple( (file.toString().split("barcode"))[1][0..1], file) }

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
    tuple val(datasetID), file("${datasetID}_filtered.fq.gz") into filtered_reads, filtered_reads_racon, filtered_reads_racon2, filtered_reads_racon3, filtered_reads_medaka, filtered_reads_nextpolish, filtered_reads_nextpolish2, filtered_reads_minimap2, filtered_reads_minimap3

    script:
        
    """
    gunzip -c ${tofilter} | NanoFilt --quality 10 --length 2500 --summary ${summary} | gzip > ${datasetID}_filtered.fq.gz
    """
}


filtered_reads.filter(){ it.get(1).countFastq() > 4}.into{ filtered_reads_canu; filtered_reads_canu_print }
//filtered_reads_canu_print.view { "name: ${it.get(1)} size ${it.get(1).countFastq()}" }

// TODO: fix path to canu
// TODO: Canu really struggles with low read counts, right now we're filtering out anything with a file size less than 100Kb
process CanuCorrect {
    publishDir "$results_path/canu"
    
    input:
    tuple val(datasetID), file(to_correct) from filtered_reads_canu
    
    output:
    tuple val(datasetID), file("${datasetID}_canu_correct/reads.correctedReads.fasta.gz") into canu_corrected_minimap, canu_corrected_convert
    
    script:

    """
    rm -rf ${datasetID}_canu_correct
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
    set str_name, file("${str_name}_racon_medaka/consensus.fasta") into racon_medaka_consensus, racon_medaka_consensus_mars
    
    script:
    str_name = tuple_pack.get(0).get(0)
    """
    medaka_consensus -i ${tuple_pack.get(0).get(1)} -d ${tuple_pack.get(1).get(1)} -o ${str_name}_racon_medaka -m $params.medaka_model

    """
}

//pilon_consensus.view { "value: $it" }
//sample_table.view { "sample_table: $it" }
read_phased_medaka_for_nextpolish = filtered_reads_nextpolish.phase(racon_medaka_consensus)


process NextPolish {
    errorStrategy 'finish'
    publishDir "$results_path/nextpolish"
   
    input:
    val tuple_pack from read_phased_medaka_for_nextpolish.filter{ it.get(1).get(1).countFasta()>=1} 

    output:
    
    set datasetID, file("${datasetID}.nextpolish.fasta") into nextpolish_consensus, nextpolish_consensus_assessment
    set datasetID, file("${datasetID}.nextpolish.fasta.stat") into nextpolish_consensus_stat, nextpolish_consensus_assessment_stat
    
    script:
    datasetID = tuple_pack.get(0).get(0)
    """
    cp ${params.nanopolish_run_config} ./run.cfg
    cp ${tuple_pack.get(1).get(1)} consensus.fasta
    echo ${tuple_pack.get(0).get(1)} > lgs.fofn
    ${params.nanopolish} run.cfg
    cp 01_rundir/genome.nextpolish.fasta ${datasetID}.nextpolish.fasta
    cp 01_rundir/genome.nextpolish.fasta.stat ${datasetID}.nextpolish.fasta.stat
    """
}


read_phased_nextpolish_for_nextpolish2 = filtered_reads_nextpolish2.phase(nextpolish_consensus)

process NextPolish2CommaThePolishing {
    errorStrategy 'finish'
    publishDir "$results_path/nextpolish2"
   
    input:
    val tuple_pack from read_phased_nextpolish_for_nextpolish2.filter{ it.get(1).get(1).countFasta()>=1} 

    output:
    
    set datasetID, file("${datasetID}.nextpolish2.fasta") into nextpolish_consensus2, nextpolish_consensus_assessment2
    set datasetID, file("${datasetID}.nextpolish2.fasta.stat") into nextpolish_consensus_stat2, nextpolish_consensus_assessment_stat2
    
    script:
    datasetID = tuple_pack.get(0).get(0)
    """
    cp ${params.nanopolish_run_config} ./run.cfg
    cp ${tuple_pack.get(1).get(1)} consensus.fasta
    echo ${tuple_pack.get(0).get(1)} > lgs.fofn
    ${params.nanopolish} run.cfg
    cp 01_rundir/genome.nextpolish.fasta ${datasetID}.nextpolish2.fasta
    cp 01_rundir/genome.nextpolish.fasta.stat ${datasetID}.nextpolish2.fasta.stat
    """
}

//pilon_consensus.view { "value: $it" }
//sample_table.view { "sample_table: $it" }
read_phased_nextpolish2 = nextpolish_consensus2.phase(sample_table)
read_phased_nextpolish2.into{nextpolish2_consensus_for_mars; nextpolish_consensus2_for_copy}

process LCPCorrection {
    errorStrategy 'finish'
    publishDir "$results_path/lcp"
   
    input:
    val tuple_pack from nextpolish2_consensus_for_mars
    //set samplePosition, sample_ID, reads1, reads2, [sangers]  from sample_table

    
    output:
    set str_name, path("${str_name}_corrected.fasta") into lcp_corrected
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    python /analysis/2021_05_06_nanopore_pipeline/PlasmidSeq/scripts/processing/suffix_array_dup_detection.py --plasmid_fasta ${tuple_pack.get(0).get(1)} --output_fasta ${str_name}_corrected.fasta
    """
}

read_phased_lcp = lcp_corrected.phase(sample_table_assessment)
read_phased_lcp.into{read_phased_lcp2}

process Rotated {
    errorStrategy 'finish'
    publishDir "$results_path/rotated"
   
    input:
    val tuple_pack from read_phased_lcp2
    
    
    output:
    set str_name, path("${str_name}_rotated.fasta"), path("${str_name}_rotated_reference.fasta")  into post_lcp_rotated
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    python /analysis/2021_05_06_nanopore_pipeline/PlasmidSeq/scripts/processing/orient_contigs.py --assembly ${tuple_pack.get(0).get(1)} --reference ${tuple_pack.get(1).get(2)} --assembly_out ${str_name}_rotated.fasta --reference_out ${str_name}_rotated_reference.fasta
    """
}


process SampleCopy {
    errorStrategy 'finish'
    publishDir "$results_path/simple_copy"
   
    input:
    val tuple_pack from nextpolish_consensus2_for_copy
    //set samplePosition, sample_ID, reads1, reads2, [sangers]  from sample_table

    
    output:
    set str_name, path("${str_name}_identity.fasta") into sample_copy
    
    script:
    str_name = tuple_pack.get(1).get(1)

    """
    cp ${tuple_pack.get(0).get(1)} ${str_name}_identity.fasta

    """
}



process AlignAndCompare {
    errorStrategy 'finish'
    publishDir "$results_path/comparison_basic"
   
    input:
    tuple sample, assembled, ref from post_lcp_rotated
    
    output:
    set sample, path("${sample}_aligned.fasta") into needleall_fasta
    
    script:

    """
    cat ${ref} ${assembled} > full.fa

    needleall -asequence ${ref} -bsequence ${assembled} -gapopen 10 -gapextend 0.5 -aformat fasta -outfile ${sample}_aligned.fasta
    

    """
}

/*
process PlasmidComparison {
    errorStrategy 'finish'
    publishDir "$results_path/comparison_pilon"
   
    input:
    val tuple_pack from read_phased_pilon2
    
    output:
    set str_name, path("${str_name}_rotated.fasta") into plasmid_comp
    
    script:
    str_name = tuple_pack.get(1).get(1)

    """

    assess_assembly.py ${tuple_pack.get(0).get(1)} ${tuple_pack.get(1).get(2)} --mode genome

    """
}*/
