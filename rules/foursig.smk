rule motif_finder:
    input:
        ref = config['ref'],
    params:
        organism = config['organism'],
        re1 = config['REEnz'][0],
        re2 = config['REEnz'][1],
        res1 = config['RESeq'][0],
        res2 = config['RESeq'][1],
        motif1 = f"scripts/fourSig/{config['organism']}_re_sites/b{len(config['RESeq'][0])}_{config['REEnz'][0]}_sites.txt",
        motif2 = f"scripts/fourSig/{config['organism']}_re_sites/b{len(config['RESeq'][1])}_{config['REEnz'][1]}_sites.txt"
    output:
        outputdir + "/fourSig/motif_finished.txt"
    conda:
        "../envs/motif.yml"
    log:
        '../logs/motiffinder.log'
    shell:
        '''
        echo "Looking for motifs ...\n"
        Rscript  ./scripts/makeReSiteFile_dyn.R {configpath}
        echo "finished motif"

         if [[ ! -f "./{params.motif1}" ]]; then
            echo "./{params.motif1} not found. Create motif. This can take some time...\n"
            if [[ {params.organism} == "hg"* ]]; then
	        perl ./scripts/fourSig/motifFinder.pl -H "{params.res2}" "{input.ref}" > "./{params.motif1}"
	    elif [[ {params.organism} == "mm10"* ]]; then
	        perl ./scripts/fourSig/motifFinder.pl -M "{params.res2}" "{input.ref}" > "./{params.motif1}"
	    else
	        "Not a valid org"  #TODO rm
	    fi
        else
            echo "Motif file exists.\n"
        fi


        if [[ ! -f "./{params.motif2}" ]]; then
           echo "./{params.motif2} not found. Create motif. This can take some time..."
           mkdir -p $(dirname {params.motif2})
           if [[ {params.organism} == "hg"* ]]; then
               (perl ./scripts/fourSig/motifFinder.pl -H "{params.res2}" "{input.ref}" > "./{params.motif2}")  
           elif [[ {params.organism} == "mm10"* ]]; then
               (perl ./scripts/fourSig/motifFinder.pl -M "{params.res2}" "{input.ref}" > "./{params.motif2}") 
           else
               echo "Not a valid org"
           fi
        else
            echo "Motif file exists."
        fi

        touch {output}
        '''
 
 
 
rule bam_to_re_tab:
    input:
        bam = rules.sort_bam.output,
        motif = rules.motif_finder.output
    params:
        organism = config['organism'],
        motif1 = f"scripts/fourSig/{config['organism']}_re_sites/b{len(config['RESeq'][0])}_{config['REEnz'][0]}_sites.txt",
        motif2 = f"scripts/fourSig/{config['organism']}_re_sites/b{len(config['RESeq'][1])}_{config['REEnz'][1]}_sites.txt"
    output:
        outputdir + '/fourSig/{sample}.tab'
    conda:
        "../envs/foursig.yml"
    log:
        '../logs/{sample}_bamtoretab.log'
    shell:
        '''
        echo "\n>>>>>>>>>>     Creating tab files (fourSig)         <<<<<<<<<<<\n"
        if [[ "{params.organism}" == "hg"* ]]; then
            echo "{params.organism} human"
            (
               perl ./scripts/fourSig/bamToReTab.pl -H -o "Y" 18 700 20 "{input.bam}" "./{params.motif1}" "./{params.motif2}" NONE > "{output}"
            ) || touch {output}
            if [ $? -ne 0 ]; then
                echo "Error creating fourSig tab file for human" >&2
            fi
        elif [[ "{params.organism}" == "mm"* ]]; then
            echo "{params.organism} m"
            (
               perl ./scripts/fourSig/bamToReTab.pl -M -o "Y" 18 700 20 "{input.bam}" "./{params.motif1}" "./{params.motif2}" NONE > "{output}"
            ) || touch {output}
            if [ $? -ne 0 ]; then
                echo "Error creating fourSig tab file for mouse" >&2
            fi
        else
            echo "No valid organism specified"
            exit 1
        fi
        if [ ! -f "{output}" ]; then
            echo "Output file not generated, exiting..." >&2
            exit 1
        fi
        '''


rule foursig_it:
    input:
        tab = rules.bam_to_re_tab.output
    output:
         outputdir + '/fourSig/{sample}_fourSigit_finished.txt'
    conda: 
        "../envs/foursig.yml"
    log:
        '../logs/{sample}_foursig.log'
    shell:
        '''
        Rscript ./scripts/run_fourSig_it.R ./{input.tab} {configpath}
        touch {output}
        '''
