rule peakC_FDR:
    input:
        wig = expand("results/{author}/basic4cseq/wig/{cond}_{rep}.wig", cond=conditionsclean, rep=config['conditionRep'], author=config['author'])
    output:
        dummy = outputdir + '/peakC/peakC_FDR_done.txt'
    conda: 
        "../envs/peakc.yml"
    shell:
        '''
        echo "\n>>>>>>>>>>     Starting peakC          <<<<<<<<<<<\n"
        Rscript ./scripts/run_peakC_FDR.R {configpath} 
        touch {output.dummy}
        '''


