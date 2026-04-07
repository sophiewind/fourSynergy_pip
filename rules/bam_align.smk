rule bwamem_align:
    input:
        sample =  inputdir + "/{sample}.fastq" ,
        ref = config['ref']
    output:
        temp(outputdir + '/alignment/{sample}_unsorted.bam')
    log:
        '../logs/{sample}_align.log'
    conda:
        "../envs/align.yml"
    shell:
        '''
        echo "\n>>>>>>>>>>     Starting Alignment          <<<<<<<<<<<\n"
        which samtools
        bwa mem -t 12 {input.ref} {input.sample} 2> {log} | samtools view -b - > {output}
        '''

rule sort_bam:
    input:
        bam = rules.bwamem_align.output
    output:
        outputdir + '/alignment/{sample}_sorted.bam'
    conda:
        '../envs/align.yml'
    log:
        '../logs/{sample}_sorting.log'
    shell:
        '''
        echo "\n>>>>>>>>>>     Sort bams          <<<<<<<<<<<\n"
        samtools sort -o {output} -@ {threads} {input.bam}
        '''

rule index_bam:
    input:
        bam = rules.sort_bam.output
    output:
        outputdir + '/alignment/{sample}_sorted.bam.bai'
    conda:
        "../envs/align.yml"
    log:
        '../logs/{sample}_indexing.log'
    shell:
        '''
        echo "\n>>>>>>>>>>     Index bams          <<<<<<<<<<<\n"
        samtools index {input} {output}
        '''
