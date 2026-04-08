rule postprocessing:
    input:
        outputdir + '/r4cker/r4cker_done.txt',
        expand('results/{author}/fourSig/{cond}_{rep}_fourSigit_finished.txt', cond=conditionsclean, rep=config['conditionRep'], author=config['author']),
        outputdir + '/r4cker/r4cker_done.txt',
        outputdir + '/peakC/peakC_FDR_done.txt',
        outputdir + '/r3cseq/r3cseq_w_done.txt',

    output:
        outputdir + "/fourSynergy_postprocessing.html",
        outputdir + "/nearbait_area.bed",
        output_dummy = outputdir + '/sia/finished.txt'
    conda:
        "../envs/ensemble.yml"
    log:
        '../logs/postprocessing.log'
    shell:
        '''
        echo "{outputdir}"
        Rscript --vanilla -e \
        'rmarkdown::render("./scripts/fourSynergy_postprocessing.Rmd", output_dir="{outputdir}",\
        params = list(config="{configpath}"))'
        touch {output.output_dummy}
        '''
