from snakemake.utils import min_version

min_version("7.0")

include_shiny = '--no-shiny' not in sys.argv

inputdir = config['inputDir']
configpath = inputdir + "/info.yaml"
outputdir = "results/" + config['author']

diff = False

conditions = [config['condition'], config['control']]


conditionsclean = [x for x in conditions if x and x.strip()]

REPS = list(map(str, config["conditionRep"]))  # macht sicher Strings
SAMPLES = [f"{c}_{r}" for c in conditionsclean for r in REPS]

print(SAMPLES)



rule all:
    input:
        outputdir + '/r4cker/r4cker_done.txt',
        outputdir + '/r3cseq/r3cseq_w_done.txt',
        outputdir + '/peakC/peakC_FDR_done.txt',
        outputdir + "/multiqc_report.html",
        expand('results/{author}/fourSig/{cond}_{rep}_fourSigit_finished.txt', cond=conditionsclean, rep=config['conditionRep'], author=config['author']),
        outputdir + "/fourSynergy_postprocessing.html",
        outputdir + "/shiny_in/.done"


include: "rules/bam_align.smk"
include: "rules/basic4cseq.smk"
include: "rules/peakC.smk"
include: "rules/r3cseq.smk"
include: "rules/r4cker.smk"
include: "rules/foursig.smk"
include: "rules/quality_control.smk"
include: "rules/postprocessing.smk"
include: "rules/copy_shiny.smk"