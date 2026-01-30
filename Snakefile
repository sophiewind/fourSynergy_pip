from snakemake.utils import min_version
min_version("7.0")

include_shiny = '--no-shiny' not in sys.argv

inputdir = config['inputDir']
configpath = inputdir + "/info.yaml"
outputdir = "results/" + config['author']

diff = False

conditions = [config['condition'], config['control']]
conditionsclean = [x for x in conditions if x and x.strip()]
print(conditionsclean)

rule all:
    input:
        outputdir + '/r4cker/r4cker_done.txt',
        outputdir + '/r3cseq/r3cseq_w_done.txt',
        outputdir + '/peakC/peakC_FDR_done.txt',
        outputdir + "/multiqc_report.html",
        outputdir + "/fourSynergy_postprocessing.html"




include: "rules/bam_align.smk"
include: "rules/basic4cseq.smk"
include: "rules/peakC.smk"
include: "rules/r3cseq.smk"
include: "rules/r4cker.smk"
include: "rules/foursig.smk"
include: "rules/quality_control.smk"
include: "rules/postprocessing.smk"
include: "rules/copy_shiny.smk"
