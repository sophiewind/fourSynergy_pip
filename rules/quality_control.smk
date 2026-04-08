rule qc_fastqc:
    input:
        fastq = expand('Datasets/{author}/{cond}_{rep}.fastq', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author']),
    output:
        zip = expand('results/{author}/qc/{cond}_{rep}_fastqc.zip', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author']),
    params:
        outputDir = outputdir + '/qc'
    log:
        '../logs/fastqc.log'
    conda:
        "../envs/fastqc.yml"
    shell:
        '''
        mkdir -p {params.outputDir}
        fastqc {input.fastq} -o {params.outputDir}
        '''

rule flagstat:
    input:
        bam = rules.sort_bam.output
    output:
        outputdir + '/alignment/{sample}_sorted_stats.txt'
    conda:
        "../envs/align.yml"
    log:
        '../logs/{sample}_flagstat.log'
    shell:
        '''
         echo "\n>>>>>>>>>>     flagstat          <<<<<<<<<<<\n"
        samtools flagstat {input} > {output}
        '''


rule multiqc:
    input:
        bam = expand('results/{author}/alignment/{cond}_{rep}_sorted.bam', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author']),
        bai = expand('results/{author}/alignment/{cond}_{rep}_sorted.bam.bai', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author']),
        qc = rules.qc_fastqc.output,
        flagstat = expand('results/{author}/alignment/{cond}_{rep}_sorted_stats.txt', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author']),
        #fastqc_report ='multiqc_fastqc.txt',
        #sequence_quality = 'fastqc_per_sequence_quality_scores_plot.txt'
    params:
        outputDir = outputdir
    log:
        '../logs/multiqc.log'
    output:
        report = outputdir + '/multiqc_report.html',
        out = outputdir + '/multiqc_data/samtools-flagstat-pct-table.txt',
        out1 = outputdir + '/multiqc_data/fastqc_per_sequence_quality_scores_plot.txt',
        out2 = outputdir + '/multiqc_data/fastqc_sequence_counts_plot.txt',
        out3 = outputdir + '/multiqc_data/multiqc_fastqc.txt'
    conda:
        '../envs/multiqc.yml'
    shell:
        '''
        multiqc {params.outputDir} -o {params.outputDir} --force
        ''' 


