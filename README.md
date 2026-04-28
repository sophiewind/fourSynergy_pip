
# fourSynergy_pip: Ensemble-based 4C-seq analysis preprocessing pipeline

Circular Chromosome Conformation Capture Sequencing (4C-seq) is a sequencing technique enabling the identification of chromatin interactions. Existing 4C-seq analysis tools are:

-   r3Cseq [@Thongjuea2013]

-   fourSig [@Williams2014]

-   peakC [@Geeven2018]

-   R.4Cker [@Raviram2016]

Our ensemble algorithm *fourSynergy* combines the strengths of all these tools to achieve superior predictive performance for interaction calling. *fourSynergy* consists of a preprocessing pipeline (*fourSynergy_pip*), an R/Bioconductor package (*fourSynergy*), and a Shiny app.

# Part I - pipeline

To start an analysis, you first need to run the pipeline. The pipeline reprocesses the 4C-seq data and runs all base tools used in our ensemble algorithm.

The pipeline is written in Snakemake and includes quality control, interaction calling, and post-processing steps to prepare the outputs of the base tools for the ensemble algorithm. The ensemble-based 4C-seq analysis can then be performed using the R package *fourSynergy*:\

<https://bioconductor.org/packages/fourSynergy>.

## Requirements

You need:

-   An indexed reference genome (BWA): <https://bio-bwa.sourceforge.net/bwa.shtml>
-   A config file (`info.yaml`)
-   FASTQ files from a 4C-seq experiment with at least two replicates
-   A conda or docker installation

The pipeline can be run via conda or docker since we offer a dockerimage.

First you need to clone the foursynergy pipeline git repository onto your computer by open a linux terminal and run the following code:

```         
 git clone https://github.com/sophiewind/fourSynergy_pip.git
```

Once you done that, change you directory:

```         
cd fourSynergy_pip/
```

## Setting up the config file

If you want to run the pipeline on test data, please download our simulated demo dataset:

```         
bash ./scripts/zenodo_download.sh
```

If you want to run the pipeline on your own data, you can use our Shiny app to generate an `info.yaml` [file](%5Bfile:\%5D(file:)%7B.uri%7D){.uri}:

<https://gerohei.shinyapps.io/config_app/>

Move the resulting `info.yaml` into the `Datasets` folder.

An exemplary config file can be found here: https://github.com/sophiewind/fourSynergy/blob/main/inst/extdata/Datasets/Demo/info.yaml

## Getting started - docker mode

Then build the docker image:

```         
docker build -t swind.foursynergy_pip .
```

Before starting the docker you need to copy your data (fastq files) and info.yaml into the foursynergy_pip folder. Therefore create a folder called "Datasets" and store your data there. If you want to run the pipeline on test data you can skip this step.

### Run on test data

If you need test data you can download a simulated dataset via Zenodo. Than just follow the tutorial by running the pipeline.

Start the pipeline with the following command. Adjust:

-   the reference genome path (before `:/ref`)

-   the dataset path / config file path (`info.yaml`) (on test data this should be fine)

The reference genome must be indexed (<https://bio-bwa.sourceforge.net/bwa.shtml>).

```         
docker run --rm \
  -v "$(pwd)":/workflow \
  -v /media/home2/share/Genomes/Mus_musculus.GRCm38/bwa_0.7.10/:/ref \  # adjust
  -w /workflow \
  swind.foursynergy_pip:latest \
  snakemake -s Snakefile --cores 30 --configfile ./Datasets/m4+4_DC1/info.yaml --use-conda --resources mem_mb=16000 -j 4
```

You can adjust resource usage and many more via the Snakemake CLI: <https://snakemake.readthedocs.io/en/stable/executing/cli.html>

If you want to continue with standard R analysis jump to "Part IIa" if you want to use the Shiny app jump to "Part IIb". 

## Getting started - conda mode

Create a snkameka conda env and activate it:

```         
conda create -n snakemake_env -c conda-forge -c bioconda snakemake
conda activate snakemake_env
```

Clone the repository:

```         
 git clone https://github.com/sophiewind/fourSynergy_pip.git
```

Change the directory:

```         
cd fourSynergy_pip/
```

Check your `info.yaml`. Paths in conda mode may differ from paths in docker mode.

Start the workflow:

```         
snakemake  -s Snakefile --cores 30 --configfile ./Datasets/m4+4_DC1/info.yaml --use-conda --resources mem_mb=16000 -j 4
```

You can add the `-n` flag for a dry run to check whether all files are created correctly.

If you want to continue with standard R analysis jump to "Part IIa" if you want to use the Shiny app jump to "Part IIb". 


# Part IIa - R analysis

# Requirements

To run *fourSynergy*, you need R Version 4.5.1 or higher.

## Installation

fourSynergy is available via Bioconductor:

```         
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("fourSynergy")
```

The latest development version from GitHub can be installed via:

```         
if (!requireNamespace("devtools", quietly=TRUE)) install.packages("devtools")
devtools::install_github("sophiewind/fourSynergy")
```

## Analyses

fourSyngery offers a variety of analyses:

-   Base tool interaction results
-   Visualization of base tool interactions
-   Ensemble interaction calling
-   Visualization of ensemble interactions
-   Differential interaction calling
-   Visualization of differential interactions
-   Highlighting of genes of interest (Note: For easier analyzation of the results, we provided a shiny app, to compare and evaluate the output.) \# Detailed documentation

# Part IIb - R shiny analysis

To analyze the data via our R Shiny UI you just need to start the Shiny App `app.R`.
A user guide can be found in this repo: `user_guide_shiny.zip`.

### docker mode 
Build the Shiny app docker image:
```
docker build -t swind.foursynergy_shiny -f ./Dockerfile.shiny .
```

Run the Shiny app:
```
docker run --rm   -v "$(pwd)":/workflow -w /workflow -p 8000:8000 swind.foursynergy_shiny
```

Access the Shiny app via you browser: `http://localhost:8000/`

### conda mode



## Documentation

Detailed documentation is available in the vignette of the R Bioconductor package: <https://bioconductor.org/packages/fourSynergy>.

For questions, feature requests, or bug reports, please contact Sophie Wind (sophie.wind\@uni-muenster.de).
