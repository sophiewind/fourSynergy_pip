rule r3cseq_window:
    input:
       bam = expand('results/{author}/alignment/{cond}_{rep}_sorted.bam', 
                     cond=conditionsclean, 
                     rep=config['conditionRep'], 
                     author=config['author'])
    output:
        outputdir + '/r3cseq/r3cseq_w_done.txt'
    conda: 
        "../envs/r3cseq.yml"
    log:
        '../logs/r3cseq.log'
    shell:
        '''
        bampath=$(dirname {input.bam[0]})
        for window_size in 2000 5000 10000; do
            (
                Rscript ./scripts/run_r3cseq_window.R ./$bampath {configpath} $window_size
            ) || echo "Error running R script with window size $window_size"
        done
        touch {output}
        '''
