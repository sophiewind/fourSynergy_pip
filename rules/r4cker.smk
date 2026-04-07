R4_BEDGRAPHS = expand(
    "results/{author}/alignment/{cond}_{rep}_sorted.bedGraph",
    cond=conditionsclean,
    rep=config["conditionRep"],
    author=config["author"]
)

R4_BEDGRAPHS_RM = expand(
    "results/{author}/alignment/{cond}_{rep}_sorted_rm_self_und.bedGraph",
    cond=conditionsclean,
    rep=config["conditionRep"],
    author=config["author"]
)

rule r4cker:
    input:
        bams = expand(
            "results/{author}/alignment/{cond}_{rep}_sorted.bam",
            cond=conditionsclean, rep=config["conditionRep"], author=config["author"]
        )
    output:
        bedgraphs = R4_BEDGRAPHS,
        bedgraphs_rm = R4_BEDGRAPHS_RM,
        done = outputdir + "/r4cker/r4cker_done.txt"
    params:
        author = config["author"],
        vpchr = config["VPchr"],
        vppos_start = config["VPpos"] - 1000,
        vppos_end   = config["VPpos"] + 1000
    conda:
        "../envs/r4cker.yml"
    shell:
        r"""
        mkdir -p results/{params.author}/r4cker

        echo -e "chr{params.vpchr}\t{params.vppos_start}\t{params.vppos_end}" \
          > results/{params.author}/r4cker/viewpoint_mask.bed

        # Create bedGraphs (names will match output list because input bam names match)
        for file in {input.bams}; do
            bedtools genomecov -bg -ibam "$file" > "${{file%.bam}}.bedGraph"
        done

        # Subtract mask -> _rm_self_und.bedGraph
        for bedgraph_file in results/{params.author}/alignment/*.bedGraph; do
            if [[ "$bedgraph_file" != *"_rm_self_und.bedGraph" ]]; then
                bedtools subtract -a "$bedgraph_file" \
                  -b results/{params.author}/r4cker/viewpoint_mask.bed \
                  > "${{bedgraph_file%.bedGraph}}_rm_self_und.bedGraph"
            fi
        done

        (Rscript ./scripts/run_r4cker.R {configpath}) || echo "R script failed"

        touch {output.done}
        """
