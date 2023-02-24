genome = config["genome"]
rule read_distro:
    input:
        "{sample}_Aligned.sortedByCoord.out_mkdup_filtered.bam"
    output:
        "{sample}_Aligned.sortedByCoord.out_mkdup_filtered_distro.txt"
    run:
        "read_distribution.py -i {input} -r ~/genomes/{genome}/rseqc/*.bed > {output}"