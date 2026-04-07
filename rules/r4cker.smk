rule r4cker:
    input:
        bams = expand('results/{author}/alignment/{cond}_{rep}_sorted.bam', cond=conditionsclean, rep=config['conditionRep'], author=config['author']),
        viewpoint = inputdir + '/viewpoint.bed'
    output:
        #outputdir + '/r4cker/out/{sample}_nearbait_highinter.bed',
        #outputdir + '/alignment/{sample}_sorted.bedGraph'
        outputdir + '/r4cker/r4cker_done.txt'
    params:
        author = config['author'],
        vpchr = config['VPchr'],
        vppos_start = config['VPpos'] - 1000,
        vppos_end = config['VPpos'] + 1000
    conda:
        "../envs/r4cker.yml"
    shell:
        '''
        echo "\n>>>>>>>>>>     Run   r4cker       <<<<<<<<<<<\n"
        echo "chr{params.vpchr}\t{params.vppos_start}\t{params.vppos_end}" > ./results/{params.author}/r4cker/viewpoint_mask.bed
        
         # Create bedGraphs
        for file in {input.bams}; do
            bedtools genomecov -bg -ibam $file > "${{file%.bam}}.bedGraph"
        done

        # Subtract viewpoint_mask.bed from bedGraphs
        for bedgraph_file in ./results/{params.author}/alignment/*bedGraph; do
        if [[ "$bedgraph_file" != *"_rm_self_und.bedGraph" ]]; then
            bedtools subtract -a $bedgraph_file -b ./results/{params.author}/r4cker/viewpoint_mask.bed > "${{bedgraph_file%.bedGraph}}_rm_self_und.bedGraph"
        fi
        done

        (
            Rscript ./scripts/run_r4cker.R {configpath}
        ) || echo "R script failed, but continuing to create r4cker_done.txt"

        touch {output}
        '''
