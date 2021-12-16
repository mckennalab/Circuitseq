#!/usr/bin/env nextflow

/*
========================================================================================
                         mckenna_lab/plasmid_seq
========================================================================================
 Process sequenced plasmids on Oxford Nanopore into assembled maps of each sequence
----------------------------------------------------------------------------------------
*/

input_fastq5_path  = Channel.fromPath(params.fast5).first()
results_path       = "results"
mix_in_rate = 0.19845440202552694
/*
 * Read in the sample table and convert entries into a channel of sample information 
 */

Channel.fromPath( params.samplesheet )
        .splitCsv(header: true, sep: '\t')
        .map{ tuple(it.position.padLeft(2,'0'), it.sample, file(it.reference), it.sangers ) }
        .into{sample_table; sample_table_assessment; sample_table_assessment2; sample_table_pre_merge;  sample_table_pre_lf; sample_table_pre_polish; sample_table_methylation}

// check if they've asked for methylation calling, if not, set it to false
if (!binding.hasVariable('params.methylation_calling')) {
    params.methylation_calling = false
}
log.info "Methylation calling: " + params.methylation_calling

// check if they've asked for QC output, if not, set it to false
if (!binding.hasVariable('params.quality_control_processes')) {
    params.quality_control_processes = false
}
log.info "Quality control output: " + params.quality_control_processes

// check if they've asked to annotate the resulting plasmid maps
if (!binding.hasVariable('params.annotate_plasmid')) {
    params.annotate_plasmid = false
}
log.info "Annotate plasmid output: " + params.annotate_plasmid
println "Project : $workflow.projectDir"
println "Git info: $workflow.repository - $workflow.revision [$workflow.commitId]"
println "Cmd line: $workflow.commandLine"
println "Manifest's pipeline version: $workflow.manifest.version"


/*
 * Initial basecalling of all reads using Guppy -- this is the most computationally costly step
 */
process GuppyBaseCalling {
    label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/guppy"

    input:
    path input_path from input_fastq5_path

    output:
    path "basecalling/pass/" into basecalled_directory
    path "basecalling/sequencing_summary.txt" into basecalling_summary_file, basecalling_summary_for_pyco   // the fastq output file path

    when:
    !params.basecalling_dir

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
    chmod -R 777 ./
    """
}

/*
 * split samples by their tagmentation barcode
 */
process GuppyDemultiplex {
    label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/guppy_demultiplex"

    input:
    path basecalled from basecalled_directory

    output:
    path "saved_data/barcoding_summary.txt" into barcoding_split_summary   // the fastq output file path
    path "saved_data/barcode**/**.fastq.gz" into fastq_gz_split_files
    
    when:
    !params.basecalling_dir
    
    script:
        
    """
    guppy_barcoder \\
        --input_path ${basecalled} \\
        --save_path saved_data \\
        --data_path ${params.barcodes} \\
        --barcode_kits MY-CUSTOM-BARCODES \\
        --front_window_size 120 \\
        --min_score_mask 30 \\
        -x $params.gpu_slot \\
        --min_score $params.barcode_min_score \\
        -q 1000000 \\
        --compress_fastq
    chmod -R o+rw ./
    """
}

/*
 * split samples by their tagmentation barcode
 */
process GuppyDemultiplexExisting {
    label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
    beforeScript 'chmod o+rw .'
    publishDir "$results_path/guppy_demultiplex"

    input:
    path basecalled from params.basecalling_dir

    output:
    path "saved_data/barcoding_summary.txt" into barcoding_split_summary_existing   // the fastq output file path
    path "saved_data/barcode**/**.fastq.gz" into fastq_gz_split_files_existing
    
    when:
    params.basecalling_dir
    
    script:
        
    """
    guppy_barcoder \\
        --input_path ${basecalled} \\
        --save_path saved_data \\
        --data_path ${params.barcodes} \\
        --barcode_kits MY-CUSTOM-BARCODES \\
        --front_window_size 120 \\
        --min_score_mask 30 \\
        -x $params.gpu_slot \\
        --min_score $params.barcode_min_score \\
        -q 1000000 \\
        --compress_fastq
    chmod -R o+rw ./
    """
}

/*
 * Run the pyco quality control step on the whole collection of reads from guppy
 */
process pycoQC {
    publishDir "$results_path/pycoQC"
    beforeScript 'chmod o+rw .'

    input:
    path basecalling_summary from basecalling_summary_for_pyco
    path barcode_summary from barcoding_split_summary

    output:
    path "pycoQC.html" into pycoQC_HTML
    
    when:
    !params.basecalling_dir
    
    script:
        
    """
    bash && \\
    python --version > python_version.txt
    pycoQC \\
        --summary_file $basecalling_summary \\
        --barcode_file $barcode_summary \\
        --html_outfile pycoQC.html
    """
}

/*
 * Take the channel of sample-split reads and send them into multiple downstream channels
 */
fastq_gz_split_files.flatten().filter(){ it.countFastq() > 1 && !(it.getParent().contains("unclassified"))}.mix(fastq_gz_split_files_existing.flatten().filter(){ it.countFastq() > 1 && !(it.getParent().contains("unclassified"))})into{ 
    guppy_demulti; 
    methylation_reads; 
    guppy_alignment; 
}

// extract the well number as an ID -- this must current result in a 2 digit code from 00 to 99
raw_reads_for_alignment = guppy_alignment.map { file -> tuple( (file.toString().split("barcode"))[1][0..1], file) }.phase(sample_table_pre_merge)

// align the reads to the known reference for downstream QC 
process AlignReadsPre {
    publishDir "$results_path/minimap_initial"
    beforeScript 'chmod o+rw .'

    input:
    tuple reads,reference from raw_reads_for_alignment // reads and reference 

    when:
    params.quality_control_processes

    output:
    
    tuple val(str_name), path("${str_name}_aligned.bam") into minimap_pre_reads_unsorted
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam") into minimap_pre_reads
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam.bai") into minimap_pre_reads_bai
    path("${str_name}.fasta.dict") into sequence_dict_pre
    path("${str_name}.fasta") into sequence_fasta_pre

    script:
    str_name = reads.get(0)

    """
    grep -v ">" ${reference.get(2)} > reference_no_header.fa
    cat ${reference.get(2)} reference_no_header.fa > duplicate_ref.fa
    cp ${reference.get(2)} ${str_name}.fasta
    samtools dict ${str_name}.fasta > ${str_name}.fasta.dict
    samtools dict duplicate_ref.fa > duplicate_ref.fa.dict
    minimap2 -ax map-ont duplicate_ref.fa ${reads.get(1)} > ${str_name}_aligned.bam
    samtools sort -o ${str_name}_sorted_mapped_reads.bam ${str_name}_aligned.bam
    samtools index ${str_name}_sorted_mapped_reads.bam
    """
}


/*
 * Chop off adapter sequences
 */
process Porechop {
    publishDir "$results_path/porechop"
    beforeScript 'chmod o+rw .'

    input:
    tuple datasetID, path(fastq) from guppy_demulti.map { file -> tuple( (file.toString().split("barcode"))[1][0..1], file) }
    output:
    tuple datasetID,"${datasetID}_porechop.fq.gz" into porechop_output, porechop_output_for_minimap, length_for_alignment
    script:
        
    """
    porechop -i ${fastq} -o ${datasetID}_porechop.fq.gz \
        --format fastq.gz \
        --end_threshold  50 --extra_end_trim 50 \
        --discard_middle --middle_threshold 80
    """
}
raw_reads_for_lf_alignment = length_for_alignment.phase(sample_table_pre_lf)



process AlignReadsPostLengthFilter {
    publishDir "$results_path/minimap_length_filter"
        beforeScript 'chmod o+rw .'

    when:
    params.quality_control_processes

    input:
        tuple reads,reference from raw_reads_for_lf_alignment // reads and reference

    output:
        tuple val(str_name), path("${str_name}_sorted_post_lf_reads.bam") into minimap_post_lf_reads
	    tuple val(str_name), path("${str_name}_sorted_post_lf_reads.bam.bai") into minimap_post_lf_reads_bai
	        path("${str_name}.fasta.dict") into sequence_dict_post_lf
		    path("${str_name}.fasta") into sequence_fasta_post_lf

    script:
        str_name = reads.get(0)

    """
    cp ${reference.get(2)} ${str_name}.fasta
    samtools dict ${str_name}.fasta > ${str_name}.fasta.dict
    minimap2 -ax map-ont ${reference.get(2)} ${reads.get(1)} | samtools sort -o ${str_name}_sorted_post_lf_reads.bam
    samtools index ${str_name}_sorted_post_lf_reads.bam
    """
    }


/*
 * Filter reads by quality
 */
 process FilterReads {
    publishDir "$results_path/filter_reads"
    beforeScript 'chmod o+rw .'

    input:
    tuple val(datasetID), file(tofilter) from porechop_output
    path summary from basecalling_summary_file
    
    output:
    tuple val(datasetID), file("${datasetID}_filtered.fq.gz") into filtered_reads, filtered_reads_racon, filtered_reads_racon2, filtered_reads_racon3, filtered_reads_medaka, filtered_reads_nextpolish, filtered_reads_nextpolish2, filtered_reads_minimap

    script:
        
    """
    gunzip -c ${tofilter} | NanoFilt --quality 10 --length 500 --summary ${summary} | gzip > ${datasetID}_filtered.fq.gz
    """
}


filtered_reads.filter(){ it.get(1).countFastq() > 4}.into{ filtered_reads_canu; filtered_reads_canu_print }

/*
 * We use Canu to make high-quality concensus sequences from the filtered reads
*/
process CanuCorrect {
    publishDir "$results_path/canu"
    beforeScript 'chmod o+rw .'

    input:
    tuple val(datasetID), file(to_correct) from filtered_reads_canu
    
    output:
    tuple val(datasetID), file("${datasetID}_canu_correct/reads.correctedReads.fasta.gz") into canu_corrected_miniasm, canu_corrected_flye, canu_corrected_minimap, canu_corrected_convert
    
    script:

    """
    rm -rf ${datasetID}_canu_correct
    mkdir ${datasetID}_canu_correct
    canu -correct \
     -p reads -d ${datasetID}_canu_correct \
     genomeSize=8k \
     corOutCoverage=200 \
     stopOnLowCoverage=2 minInputCoverage=2 \
     -nanopore ${to_correct} 
    """
}

/*
 * Assemble corrected reads with Flye
 */
process Flye {
    publishDir "$results_path/flye"
    beforeScript 'chmod o+rw .'

    input:
    tuple val(datasetID), path(corrected_reads) from canu_corrected_flye.filter(){it.get(1).countLines() >= 20}

    output:
    tuple val(datasetID), path("${datasetID}_assembly/assembly.fasta") into flye_assembly
    
    script:
    """
    flye -g 10k --nano-corr \
         ${corrected_reads} \
         -o ${datasetID}_assembly \
         --min-overlap 1000

    """
}

/*
 * Assemble corrected reads with Miniasm 
 */
process Miniasm {
    publishDir "$results_path/miniasm"
    beforeScript 'chmod o+rw .'

    input:
    tuple val(datasetID), path(corrected_reads) from canu_corrected_miniasm

    output:
    tuple val(datasetID), path("${datasetID}_overlaps.gfa") into miniasm_overlaps
    
    script:
    """
    minimap2 -x ava-ont \
         ${corrected_reads} \
         ${corrected_reads} > ${datasetID}_overlaps.paf

    miniasm \
    -o 500 -I 0.9 -F 0.9 \
    -f ${corrected_reads} ${datasetID}_overlaps.paf > ${datasetID}_overlaps.gfa
    """
}

/*
 * We need to 'phase', or line-up, the corrected reads with their assembly results
 */
read_phased = canu_corrected_convert.phase(miniasm_overlaps)

/*
 * Check for any really poor quality assemblies
 */
process ConvertGraph {
    errorStrategy 'finish'
    publishDir "$results_path/convert_graph"
    beforeScript 'chmod o+rw .'

    input:
    val phased from read_phased.filter(){it.get(1).get(1).countLines() > 1} 
    
    output:
    tuple val(str_name), path("${str_name}_contigs.fasta") into miniasm_assembly

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


/*
 * We need to 'phase', or line-up, the corrected reads with their assembly results
 */
assembly_phase = flye_assembly.phase(miniasm_assembly, remainder: true)

/*
 * Check for any really poor quality assemblies
 */
process AssessAssemblyApproach {
    errorStrategy 'finish'
    publishDir "$results_path/assembly_choice"
    beforeScript 'chmod o+rw .'

    input:
    val phased_assemblies from assembly_phase
    
    output:
    tuple val(str_name), path("${str_name}_${method}_assembly.fasta") into fasta_graph

    script:

    str_name = phased_assemblies.get(0).get(0)
    myFlye = (phased_assemblies.get(0) == null) ? 'mydefaultvalue' : phased_assemblies.get(0).get(1)
    myMiniasm = (phased_assemblies.get(1) == null) ? 'mydefaultvalue' : phased_assemblies.get(1).get(1)
    method = (phased_assemblies.get(1) == null) ? "flye" : "miniasm"
    """
    if [ -f "${myMiniasm}" ]; then
        cp ${myMiniasm} ${str_name}_${method}_assembly.fasta
    else 
        cp ${myFlye} ${str_name}_${method}_assembly.fasta
    fi
    """
}


/*
 * phase the reads to the filtered assembly from the last step
 */
read_phased_racon = filtered_reads_racon.phase(fasta_graph)

/*
 * Perform the first round of polishing with Racon
 *
 * TODO: install racon the right way
 */
process RaconPolish {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_racon.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon1.fasta") into racon_corrected
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -u -m 8 -x -6 -g -8 -w 500 -t 1 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon1.fasta
    """
}


/*
 * phase the racon corrected reference to the reads again
 */
read_phased_racon2 = filtered_reads_racon2.phase(racon_corrected)

/*
 * Racon polish again
 *
 * TODO: install racon the right way
 *
 */
process RaconPolish2 {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_racon2.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon2.fasta") into racon_corrected2
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -u -m 8 -x -6 -g -8 -w 500 -t 1 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon2.fasta
    """
}

/*
 * again, phase the racon corrected reference to the reads
 */
read_phased_racon3 = filtered_reads_racon3.phase(racon_corrected2)

/*
 * Racon polish again
 *
 * TODO: install racon the right way
 *
 */
process RaconPolish3 {
    errorStrategy 'finish'
    publishDir "$results_path/racon_polish"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_racon3.filter{ it.get(1).get(1).countFasta()>=1} 

    
    output:
    tuple str_name, path("${str_name}_racon3.fasta") into racon_corrected3
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    minimap2 -ax map-ont ${tuple_pack.get(1).get(1)} ${tuple_pack.get(0).get(1)} > ${tuple_pack.get(0).get(0)}_mapping.sam

    racon -u -m 8 -x -6 -g -8 -w 500 -t 1 ${tuple_pack.get(0).get(1)} ${tuple_pack.get(0).get(0)}_mapping.sam ${tuple_pack.get(1).get(1)} > ${tuple_pack.get(0).get(0)}_racon3.fasta
    """
}

/*
 * our last phase of the polished, racon corrected reference to the reads 
 */
read_phased_medaka = filtered_reads_medaka.phase(racon_corrected3)

/*
 * perform a Medaka Consensus of the reference
 */
process MedakaConsensus {
    errorStrategy 'finish'
    publishDir "$results_path/medaka_consensus"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_medaka.filter{ it.get(1).get(1).countFasta()>=1} 
    
    output:
    set str_name, file("${str_name}_racon_medaka/consensus.fasta") into racon_medaka_consensus, racon_medaka_consensus_mars
    
    script:
    str_name = tuple_pack.get(0).get(0)
    """
    chmod -R a+rw ./
    medaka_consensus -i ${tuple_pack.get(0).get(1)} -d ${tuple_pack.get(1).get(1)} -o ${str_name}_racon_medaka -m r941_min_high_g360
    chmod -R a+rw ./
    """
}

/*
 * phase our nextpolish reads to the medaka consensus
 */
read_phased_medaka_for_nextpolish = filtered_reads_nextpolish.phase(racon_medaka_consensus)

/*
 * develop a nextpolish consensus
 */
process NextPolish {
    errorStrategy 'finish'
    publishDir "$results_path/nextpolish"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_medaka_for_nextpolish.filter{ it.get(1).get(1).countFasta()>=1} 
    path config from params.nanopolish_run_config
    output:
    
    set datasetID, file("${datasetID}.nextpolish.fasta") into nextpolish_consensus, nextpolish_consensus_assessment
    set datasetID, file("${datasetID}.nextpolish.fasta.stat") into nextpolish_consensus_stat, nextpolish_consensus_assessment_stat
    
    script:
    datasetID = tuple_pack.get(0).get(0)

    """
    cp ${tuple_pack.get(1).get(1)} consensus.fasta
    echo ${tuple_pack.get(0).get(1)} > lgs.fofn
    ${params.nanopolish} ./run.cfg
    cp 01_rundir/genome.nextpolish.fasta ${datasetID}.nextpolish.fasta
    cp 01_rundir/genome.nextpolish.fasta.stat ${datasetID}.nextpolish.fasta.stat
    """
}

/*
 * phase reads to the polished reference
 */
read_phased_nextpolish_for_nextpolish2 = filtered_reads_nextpolish2.phase(nextpolish_consensus)

/*
 * Polish the reference again with Nextpolish
 */
process NextPolish2CommaThePolishing {
    errorStrategy 'finish'
    publishDir "$results_path/nextpolish2"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_nextpolish_for_nextpolish2.filter{ it.get(1).get(1).countFasta()>=1} 
    path config from params.nanopolish_run_config

    output:
    
    set datasetID, file("${datasetID}.nextpolish2.fasta") into nextpolish_consensus2, nextpolish_consensusTest, nextpolish_consensus_assessment2
    set datasetID, file("${datasetID}.nextpolish2.fasta.stat") into nextpolish_consensus_stat2, nextpolish_consensus_assessment_stat2
    
    script:
    datasetID = tuple_pack.get(0).get(0)
    """
    cp ${tuple_pack.get(1).get(1)} consensus.fasta
    echo ${tuple_pack.get(0).get(1)} > lgs.fofn
    ${params.nanopolish} run.cfg
    cp 01_rundir/genome.nextpolish.fasta ${datasetID}.nextpolish2.fasta
    cp 01_rundir/genome.nextpolish.fasta.stat ${datasetID}.nextpolish2.fasta.stat
    """
}

/*
 * phase the polished assembly to our sample information, and make two copies
 */
read_phased_nextpolish2 = nextpolish_consensus2.phase(sample_table)
read_phased_nextpolish2.into{nextpolish2_consensus_for_LCP; nextpolish_consensus2_for_copy}

/*
 * try to address issues with large duplicated segments in the resulting plasmids
 */
process LCPCorrection {
    errorStrategy 'finish'
    publishDir "$results_path/lcp"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from nextpolish2_consensus_for_LCP
    
    output:
    set str_name, path("${str_name}_corrected.fasta") into lcp_corrected
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    python /plasmidseq/scripts/processing/suffix_array_dup_detection.py --plasmid_fasta ${tuple_pack.get(0).get(1)} --output_fasta ${str_name}_corrected.fasta
    """
}

/*
 * phase this new reference with the sample table
 */
read_phased_lcp = lcp_corrected.phase(sample_table_assessment)

/*
 * attempt to align the known reference and our assembly for downstream comparison
 */
process Rotated {
    errorStrategy 'finish'
    publishDir "$results_path/rotated"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from read_phased_lcp
    
    output:
    set str_name, path("${str_name}_rotated.fasta"), path("${str_name}_rotated_reference.fasta")  into rotated_reference_align, rotated_reference_minimap, rotated_reference_minimap2, rotated_reference_eval, rotated_reference_methyl
    
    script:
    str_name = tuple_pack.get(0).get(0)

    """
    python /plasmidseq/scripts/processing/orient_contigs.py --assembly ${tuple_pack.get(0).get(1)} --reference ${tuple_pack.get(1).get(2)} --assembly_out ${str_name}_rotated.fasta --reference_out ${str_name}_rotated_reference.fasta
    """
}

/*
 * take the orginal barcode split reads and extract a sample name for methylation analysis
 */
methylation_reads_samples = methylation_reads.map { file -> tuple( (file.toString().split("barcode"))[1][0..1], file) }

/*
 * Given the barcode split read names, extract their signals from the Fast5 file
 */
process Fast5Subset {
    errorStrategy 'finish'
    publishDir "$results_path/methylation"
    beforeScript 'chmod o+rw .'

    when:
    params.methylation_calling

    input:
    path input_path from input_fastq5_path
    set str_name,reads from methylation_reads_samples
    
    output:
    set str_name, path("${str_name}_fast5_subset") into fast5_subset
    
    script:
    """
    zcat ${reads} | awk '{if(NR%4==1) print \$1}' | sed -e "s/^@//" | cat > read_ids.txt
    fast5_subset --input ${input_path} --save_path ${str_name}_fast5_subset --read_id_list read_ids.txt --filename_base subset --batch_size 10000 --recursive
    """
}

/*
 * Create two read piles for each of the methylation types
 */
fast5_phased_base = rotated_reference_methyl.phase(fast5_subset)
fast5_phased_base.set{fast5_phased}

/*
 * Call methylation using an older guppy model and modPhred
 */
process OGMethylationCalling {
    label (params.GPU == "ON" ? 'with_gpus': 'with_cpus')
    beforeScript 'chmod o+rw .'
    errorStrategy 'finish'
    publishDir "$results_path/methylation"
    maxForks 1

    when:
    params.methylation_calling

    input:
    val tuple_pack from fast5_phased
    
    output:
    path("modPhred/${str_name}") into five_methyl2
    script:
    str_name = tuple_pack.get(0).get(0)
    """
    /og_methyl/ont-guppy/bin/guppy_basecaller -x cuda:0 -c dna_r9.4.1_450bps_modbases_dam-dcm-cpg_hac.cfg --compress_fastq --fast5_out --disable_pings -ri ${tuple_pack.get(1).get(1)} -s basecalled_mod
    /og_methyl/modPhred/run -f ${tuple_pack.get(0).get(1)} -o modPhred/${str_name} -i basecalled_mod/workspace --minModFreq 0.0 
    """
}
/*
 * Create a reference file with the name from the sample table
 */
process ReferenceCopy {
    errorStrategy 'finish'
    publishDir "$results_path/reference_copy"
    beforeScript 'chmod o+rw .'

    input:
    val tuple_pack from nextpolish_consensus2_for_copy
    
    output:
    set str_name, path("${sample_ID}_${str_name}_identity.fasta") into sample_copy
    
    script:
    str_name = tuple_pack.get(1).get(1)
    sample_ID = tuple_pack.get(0).get(0)

    """
    cp ${tuple_pack.get(0).get(1)} ${sample_ID}_${str_name}_identity.fasta

    """
}

/*
 * Create an alignment of the reads to the final reference using minimap
 */
reference_and_reads = porechop_output_for_minimap.phase(rotated_reference_minimap)
reference_and_reads.into{reference_and_reads_align; reference_and_reads_align2}

/*
 * Align the original reads back to the resulting assembled reference 
 */

process AlignReads {
    publishDir "$results_path/minimap_final_porechop"
    beforeScript 'chmod o+rw .'

    input:
    tuple reads,reference from reference_and_reads_align // reads and reference 

    output:
    tuple val(str_name), path("${str_name}_aligned.bam") into minimap_reads_unsorted
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam") into minimap_reads
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam.bai") into minimap_reads_bai
    path("${str_name}.fasta.dict") into sequence_dict
    path("${str_name}.fasta") into sequence_fasta
    
    
    script:
    str_name = reads.get(0)

    """
    grep -v ">" ${reference.get(1)} > reference_no_header.fa
    cat ${reference.get(1)} reference_no_header.fa > duplicate_ref.fa
    cp ${reference.get(1)} ${str_name}.fasta
    samtools dict ${str_name}.fasta > ${str_name}.fasta.dict
    samtools dict duplicate_ref.fa > duplicate_ref.fa.dict
    minimap2 -ax map-ont duplicate_ref.fa ${reads.get(1)} > ${str_name}_aligned.bam
    samtools sort -o ${str_name}_sorted_mapped_reads.bam ${str_name}_aligned.bam
    samtools index ${str_name}_sorted_mapped_reads.bam
    """
}

reference_and_reads_nanofilt = filtered_reads_minimap.phase(rotated_reference_minimap2)
reference_and_reads_nanofilt.set{reference_and_reads_align_nanofilt}
/*
 * QC process to check out how reads align after the nanofilter step
*/
process AlignReadsNanofilter {
    publishDir "$results_path/minimap_final_nanofilter"
    beforeScript 'chmod o+rw .'
    
    when:
    params.quality_control_processes

    input:
    tuple reads,reference from reference_and_reads_align_nanofilt // reads and reference 

    output:
    tuple val(str_name), path("${str_name}_aligned.bam"),path("${str_name}_duplicate_ref.fa") into minimap_reads_unsorted_nanofilter
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam") into minimap_reads_nanofilter
    tuple val(str_name), path("${str_name}_sorted_mapped_reads.bam.bai") into minimap_reads_bai_nanofilter
    path("${str_name}.fasta.dict") into sequence_dict_nanofilter
    path("${str_name}.fasta") into sequence_fasta_nanofilter
    
    
    script:
    str_name = reads.get(0)

    """
    grep -v ">" ${reference.get(1)} > reference_no_header.fa
    cat ${reference.get(1)} reference_no_header.fa > ${str_name}_duplicate_ref.fa
    cp ${reference.get(1)} ${str_name}.fasta
    samtools dict ${str_name}.fasta > ${str_name}.fasta.dict
    samtools dict ${str_name}_duplicate_ref.fa > ${str_name}_duplicate_ref.fa.dict
    minimap2 -ax map-ont ${str_name}_duplicate_ref.fa ${reads.get(1)} > ${str_name}_aligned.bam
    samtools sort -o ${str_name}_sorted_mapped_reads.bam ${str_name}_aligned.bam
    samtools index ${str_name}_sorted_mapped_reads.bam
    """
}

/*
 * Run a script that checks how well our aligned BAM files do on a number of contamination metrics 
*/
process AssessContamination {
    publishDir "$results_path/contamination"
    beforeScript 'chmod o+rw .'
    
    when:
    params.quality_control_processes

    input:
    tuple sample, aligned_bam, reference from minimap_reads_unsorted_nanofilter

    output:
    tuple sample, "${sample}_all_levels.txt" into contamination_all
    path("${sample}_point_estimate.txt") into contamination

    script:

    """
    python /plasmidseq/scripts/contamination/2021_11_24_estimate_contamination.py \
    --bamfile ${aligned_bam} \
    --reference ${reference} \
    --sample ${sample} \
    --plasmid_mix_in_rate ${mix_in_rate} \
    --sample_range ${sample}_all_levels.txt \
    --sample_estimate ${sample}_point_estimate.txt
    """
}


/*
 * aggregate all the plasmid comparisons into a single file
 */
process ContaminationAggregation {
    errorStrategy 'finish'
    publishDir "$results_path/agg_contam"
    beforeScript 'chmod o+rw .'
    
    when:
    params.quality_control_processes

    input:
    file point_estimate from contamination.toList()
    
    output:

    file('all_contamination_stats.txt')
    
    script:

    """
    echo "assembly\tcontamination" > header.txt
    cat header.txt ${point_estimate.collect().join(" ")} > all_contamination_stats.txt
    """
}

/*
 * Create an alignment of the reference and the known (provided) plasmid map
 */

read_phased_for_alignment = rotated_reference_align.phase(sample_table_assessment2)

process AlignReferences {
    errorStrategy 'finish'
    publishDir "$results_path/comparison_basic"
    beforeScript 'chmod o+rw .'
    
    when:
    params.quality_control_processes
    
    input:
    val tuple_pack from read_phased_for_alignment
    
    output:
    set sample_ID, path("${sample_ID}_aligned.fasta") into needleall_fasta
    
    script:
    known_ref = tuple_pack.get(1).get(2)
    sample_ID = tuple_pack.get(0).get(0)
    assembled = tuple_pack.get(0).get(2)
    """
    cat ${known_ref} ${assembled} > full.fa

    needleall -asequence ${known_ref} -bsequence ${assembled} -gapopen 10 -gapextend 0.5 -aformat fasta -outfile ${sample_ID}_aligned.fasta
    """
}

/*
 * assess the map-reference alignment
 */
process PlasmidComparison {
    errorStrategy 'finish'
    publishDir "$results_path/plasmid_comp"
    beforeScript 'chmod o+rw .'

    when:
    params.quality_control_processes

    input:
    tuple sample, assembled, ref from rotated_reference_eval
    
    output:
    path("${sample}_rotated.stats") into plasmid_comp
    
    script:

    """
    assess_assembly.py ${assembled} ${ref} --mode replicon > ${sample}_rotated.stats
    """
}

/*
 * aggregate all the plasmid comparisons into a single file
 */
process PlasmidComparisonCollection {
    errorStrategy 'finish'
    publishDir "$results_path/plasmid_stat"
    beforeScript 'chmod o+rw .'
    
    when:
    params.quality_control_processes

    input:
    file stats from plasmid_comp.toList()
    
    output:
    file('all_rotated.stats')
    
    script:

    """
    echo "assembly\treplicon_name\tlength\tcontig_name\tcontiguity\tidentity\tmax_indel" > header.txt
    cat header.txt ${stats.collect().join(" ")} > all_rotated.stats
    """
}

/*
 * aggregate all the plasmid comparisons into a single file
 */
process AnnotatePlasmid {
    errorStrategy 'finish'
    publishDir "$results_path/annotated"
    beforeScript 'chmod o+rw .'
    
    when:
    params.annotate_plasmid

    input:
    tuple name,reference from sample_copy.filter(){ it.get(1).countFasta() == 1}
    
    output:
    path("${name}_annotations") into annotation_dir
    
    script:

    """
    mkdir ${name}_annotations 
    plannotate batch -i ${reference} -o ${name}_annotations -f ${name} -b /plasmidseq/pLannotate/BLAST_dbs/

    """
}