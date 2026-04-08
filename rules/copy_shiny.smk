rule copy_shiny:
    input:
        bams=expand("results/{author}/alignment/{cond}_{rep}_sorted.bam",
                    cond=conditionsclean, rep=config["conditionRep"], author=config["author"]),
        bais=expand("results/{author}/alignment/{cond}_{rep}_sorted.bam.bai",
                    cond=conditionsclean, rep=config["conditionRep"], author=config["author"]),
        bedgraphs=expand("results/{author}/alignment/{cond}_{rep}_sorted.bedGraph",
                         cond=conditionsclean, rep=config["conditionRep"], author=config["author"]),
        basic=expand("results/{author}/basic4cseq/stats/{cond}_{rep}_stats.txt",
                     cond=conditionsclean, rep=config["conditionRep"], author=config["author"]),
        nearbait="results/{author}/nearbait_area.bed".format(author=config["author"]),
        sia_dummy="results/{author}/sia/finished.txt".format(author=config["author"]),
        mq_files=expand([
            "results/{author}/multiqc_data/multiqc_fastqc.txt",
            "results/{author}/multiqc_data/fastqc_sequence_counts_plot.txt",
            "results/{author}/multiqc_data/fastqc_per_sequence_quality_scores_plot.txt",
            "results/{author}/multiqc_data/samtools-flagstat-pct-table.txt",
        ], author=config["author"]),
        r4 = outputdir + "/r4cker/r4cker_done.txt",
        pc_FDR = outputdir + "/peakC/peakC_FDR_done.txt",
        r3c = outputdir + "/r3cseq/r3cseq_w_done.txt",

        foursig = expand(outputdir + "/fourSig/{sample}_fourSigit_finished.txt", sample=SAMPLES)
    output:
        done = outputdir + "/shiny_in/.done"
    log:
        "../logs/copy_shiny.log"
    shell:
        r"""
        mkdir -p {outputdir}/shiny_in

        for f in {input.bams} {input.bais} {input.bedgraphs} {input.basic}; do
            cp -pf "$f" {outputdir}/shiny_in/
        done

        cp -pf {input.nearbait} {outputdir}/shiny_in/nearbait_area.bed
        cp -rf $(dirname {input.sia_dummy})/* {outputdir}/shiny_in/

        for file in {input.mq_files}; do
            cp -pf "$file" {outputdir}/shiny_in/"$(basename "$file")"
        done

        touch {output.done}
        """