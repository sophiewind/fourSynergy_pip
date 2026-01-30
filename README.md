# fourSynergy_pip
*fourSynergy* is an ensemble algorithm that leverages synergies between
the base tools *r3cseq*, *fourSig*, *peakC*, and *R.4cker*. To this end,
all of these tools need to be run, which can be achieved using the
*fourSynergy* pipeline. The pipeline, written in *snakemake*, includes
quality control, interaction calling, and post-processing steps to
prepare the output of the base tools for the ensemble algorithm. The ensemble based 4C-seq analysis can than be performed in using R package fourSynergy (https://github.com/sophiewind/fourSynergy).

## Getting started
You need:
- indexed reference genome (https://bio-bwa.sourceforge.net/bwa.shtml)
- config file
- viewpoint.bed
- FASTQ files 4C-seq experiment with at least two replicates each


## Snakemake settings
You can adjust how many resources are used by snakemake via CLI (https://snakemake.readthedocs.io/en/stable/executing/cli.html).
