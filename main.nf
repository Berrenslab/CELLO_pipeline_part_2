#!/usr/bin/env nextflow

// channels
// fastqs 
fastqs = Channel.fromPath("${launchDir}/intermediates/barcode_*.fastq")
    .filter(file -> file.size() > 50000) // filter for fastq files larger than 50KB

// internal adaptor filter 
//adaptor_filter_script = Channel.fromPath("bin/internal_adaptor_filter_dT100k.Rmd")
dt_threshold_rds = Channel.fromPath("${launchDir}/intermediates/adaptor_dT_threshold.rds")
tso_threshold_rds = Channel.fromPath("${launchDir}/intermediates/adaptor_TSO_threshold.rds")

// minimap - print path of reference
reference_index = Channel.fromPath(params.minimap_reference_index)

// grouping rmd 
//grouping_script = Channel.fromPath("bin/grouping.Rmd")
barcode_rds = Channel.fromPath("${launchDir}/intermediates/barcode_*.fastq.rds")
    .map {rds -> tuple(rds.baseName.tokenize('.')[0], rds)}

// err_corr rmd 
//errcorr_script = Channel.fromPath("bin/errorcorrect.Rmd")
barcode_rds_err = Channel.fromPath("${launchDir}/intermediates/barcode_*.fastq.rds")
    .map {rds -> tuple(rds.baseName.tokenize('.')[0], rds)}

all_fastq = Channel.fromPath("${launchDir}/output/${params.experiment_name}_less20kb.fastq")

process dT_adaptor_filter {
    clusterOptions '--job-name=dt_internal_filter'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    // input
    input:    
    path fastq_file
    path dt_threshold 

    output: 
    tuple val(fastq_file.baseName), 
    path("${fastq_file}_adaptor_dT.middle.rds")

    script:
    """
    echo $fastq_file
    echo dT
    singularity exec -B $params.path $params.singularity R --vanilla -e "
    rmarkdown::render('${baseDir}/bin/internal_adaptor_filter_dT100k.Rmd', knit_root_dir = '\$PWD', intermediates_dir = '\$PWD', 
    params = list(barcode = '${fastq_file}', adaptor.type = 'dT'), output_file = '${launchDir}/output/internal_adaptor_filter_dT_${fastq_file}.html')"
    """
}

process TSO_adaptor_filter {
    clusterOptions '--job-name=tso_internal_filter'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    // input
    input: 
    path fastq_file
    path tso_threshold 

    output: 
    tuple val(fastq_file.baseName), 
    path("${fastq_file}_adaptor_TSO.middle.rds")


    script:
    """
    echo $fastq_file
    echo TSO
    singularity exec -B $params.path $params.singularity R --vanilla -e "
    rmarkdown::render('${baseDir}/bin/internal_adaptor_filter_dT100k.Rmd', knit_root_dir = '\$PWD', intermediates_dir = '\$PWD',
    params = list(barcode = '${fastq_file}', adaptor.type = 'TSO'), output_file = '${launchDir}/output/internal_adaptor_filter_TSO_${fastq_file}.html')"
    """
}

process align {
    clusterOptions '--job-name=minimap'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    input: 
    path fastq_file 
    path genome

    output: 
    tuple val(fastq_file.baseName), 
    path("${fastq_file.baseName}.sam")

    script:
    """
    module load minimap2/2.17
    out=$fastq_file
    minimap2 -ax map-ont $genome $fastq_file > \${out%fastq}sam
    """
}

process grouping{
    clusterOptions '--job-name=grouping'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    input: 
    tuple val(id), path(dt_rds), path(tso_rds), path(sam), path(fastq_rds) 

    output: 
    path "${id}.groups.*.rds"

    script:
    """
    echo $id
    singularity exec -B $params.path $params.singularity R --vanilla -e "rmarkdown::render('${baseDir}/bin/grouping.Rmd', 
   knit_root_dir = '\$PWD' , intermediates_dir = '\$PWD', params = 
  list(barcode = '$id', fastq_rds = '$fastq_rds' , sam = '$sam', dt_middle_rds = '$dt_rds', tso_middle_rds = '$tso_rds'), output_file = '${launchDir}/output/grouping_${id}.html')"
    """
}

process err_corr{
    clusterOptions '--job-name=grouping'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    input: 
    tuple val(id), path(group_rds), path(fastq_rds), path(all)

    output:
    path "barcode_*_*_corrected_all.fastq"
//  

    script:
    """
    echo $all 
    singularity exec -B $params.path $params.singularity R --vanilla -e "rmarkdown::render('${baseDir}/bin/errorcorrect.Rmd', 
   knit_root_dir = '\$PWD' , intermediates_dir = '\$PWD', params = 
  list(barcode = '$id', group_rds = '$group_rds', fastq_rds = '$fastq_rds'), output_file = '${launchDir}/output/error_corr_${group_rds}.html')"

    """

}

process corrected_merge{
    publishDir "${launchDir}/output/", mode: 'move'
    clusterOptions '--job-name=corr_merge'
    maxRetries 2
    errorStrategy { task.attempt <= 2 ? 'retry' : 'finish' }

    input:
    file correct_fastqs

    output: 
    path "corrected_barcode_*.fastq"

    script:
    """ 
    for barcode in {1..96}; do
        cat barcode_\${barcode}_*_corrected_all.fastq > corrected_barcode_\${barcode}.fastq
    done
    
    """
}


// Define the workflow
workflow {
    log.info """\
                                                                    ▐         ▐  
       §                                                         ▞▀▘▜▀ ▞▀▖▛▀▖ ▜▀ ▌  ▌▞▀▖
     *--= :                                                      ▝▀▖▐ ▖▛▀ ▙▄▘ ▐ ▖▐▐▐ ▌ ▌ 
    -+@.@%@                                                      ▀▀  ▀ ▝▀▘▌    ▀  ▘▘ ▝▀ 
          @+@                                
            @:+.  :+*=                       
             .%=@§     §@                    
            @  § :*@     §                   ┌───────────────────────────────────────────────────────────────┐
           %     % --@   +@                  │                                                               │
           §      ==  @@  *                  │               ▀▀█    ▀▀█                                      │
           §@       =   .@@ @   §            │  ▄▄▄    ▄▄▄     █      █     ▄▄▄           ▄▄▄    ▄▄▄    ▄▄▄▄ │
            @ @@:-   @ @  @ @§@. @@          │ █▀  ▀  █▀  █    █      █    █▀ ▀█         █   ▀  █▀  █  █▀ ▀█ │
              .   .@:= = :          @*       │ █      █▀▀▀▀    █      █    █   █   ▀▀▀    ▀▀▀▄  █▀▀▀▀  █   █ │
                   =@ @   +@ -        @      │ ▀█▄▄▀  ▀█▄▄▀    ▀▄▄    ▀▄▄  ▀█▄█▀         ▀▄▄▄▀  ▀█▄▄▀  ▀█▄██ │
                    @-§    %  -       .      │                                                             █ │
                  .*@        .+ :     @      │                                                             ▀ │
                    *          % @   @+      └───────────────────────────────────────────────────────────────┘ 
                     @          *  @@:       
                     *:@        %@ .         
                       - @@@@@+  §           
                                             
    
                                            Singularity: $params.singularity
                Reference genome: $params.minimap_reference_index
                                                    Experiment: $params.experiment_name
                                                        *Transposomic slay*
  """

    // adaptor filter and mapping at once 

    id_dt_tso_sam_fastqrds = dT_adaptor_filter(fastqs, dt_threshold_rds.first()) 
    .join(TSO_adaptor_filter(fastqs, tso_threshold_rds.first()))
    .join(align(fastqs, reference_index.first()))
    .join(barcode_rds) 

    grouping(id_dt_tso_sam_fastqrds).flatten()
    .map {group -> tuple(group.baseName.tokenize('.')[0], group)}
    .combine(barcode_rds_err, by: 0 )
    .combine(all_fastq.first())
    .set {grouping_out}


    corrected_merge(err_corr(grouping_out).collect())


}


