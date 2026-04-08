rule basic4cseq:
    input:
        bam = rules.sort_bam.output
    output:
        outputdir + '/basic4cseq/wig/{sample}.wig',
        outputdir + '/basic4cseq/stats/{sample}_stats.txt'
    params:
        author = config['author'],
    log:
        '../logs/{sample}_basic4cseq.log'
    conda:
        "../envs/basic4cseq.yml"
    shell:
        '''
        echo "\n>>>>>>>>>>     Run basic4cseq          <<<<<<<<<<<\n"
        Rscript ./scripts/create_wigs.R ./{input.bam} {configpath}
        '''
        
