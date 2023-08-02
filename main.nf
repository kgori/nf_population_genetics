nextflow.enable.dsl=2

/** Input VCF file containing all data */
params.vcfFile

/** Comma-separated list of files listing samples to remove from the VCF */
params.sampleList

/** Files naming outgroup samples to remove from PCA */
params.outgroupList

params.outDir

process concat_filter_lists {
    input:
    path sample_list

    output:
    path 'filter_list.txt', emit: filter_list

    script:
    """
    cat ${sample_list} > filter_list.txt
    """
}

process concat_outgroup_lists {
    input:
    path outgroup_list

    output:
    path 'outgroup_list.txt', emit: outgroup_list

    script:
    """
    cat ${outgroup_list} > outgroup_list.txt
    """
}

process remove_bad_samples {
    cpus 4

    input:
    path vcf
    path filter_list

    output:
    path "${vcf.getBaseName(2)}_filtered.vcf.gz"

    script:
    """
    bcftools view --threads 4 -S ^"${filter_list}" --force-samples "${vcf}" -Oz -o "${vcf.getBaseName(2)}_filtered.vcf.gz"
    """
}

process select_transversion_sites {
    cpus 4

    input:
    path vcf

    output:
    tuple path("${vcf.getBaseName(2)}_transversions.vcf.gz"), path("${vcf.getBaseName(2)}_transversions.vcf.gz.csi")

    script:
    """
    extract_transversions.sh "${vcf}" "${vcf.getBaseName(2)}_transversions.vcf.gz" ${task.cpus / 2}
    bcftools index --threads ${task.cpus} "${vcf.getBaseName(2)}_transversions.vcf.gz"
    """
}

process linkage_pruning {
    cpus 4

    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_pruned.vcf.gz"), path("${vcf.getBaseName(2)}_pruned.vcf.gz.csi")

    script:
    """
    make_annotation_bed.sh "${vcf}" > annotation.bed.gz
    bcftools annotate --threads ${task.cpus} -c CHROM,FROM,TO,ID -a annotation.bed.gz -o annotated.vcf.gz "${vcf}"
    plink --threads ${task.cpus} \
      --indep-pairwise 50 5 0.5 \
      --const-fid \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf annotated.vcf.gz \
      --out plink_pruning
    plink --extract plink_pruning.prune.in \
      --const-fid \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf annotated.vcf.gz \
      --recode vcf \
      --out "${vcf.getBaseName(2)}_pruned"
    bcftools query -l "${vcf.getBaseName(2)}_pruned.vcf" | sed 's/0_//' > reheader.txt &&
      bcftools reheader --threads ${task.cpus} -s reheader.txt -o tmp "${vcf.getBaseName(2)}_pruned.vcf" &&
      bcftools view --threads ${task.cpus} -o "${vcf.getBaseName(2)}_pruned.vcf.gz" -Oz tmp && rm tmp
    bcftools index --threads ${task.cpus} "${vcf.getBaseName(2)}_pruned.vcf.gz"
    rm annotation.bed.gz
    rm plink_pruning*
    rm reheader.txt
    """
}

process remove_outgroups {
    cpus 4

    input:
    tuple path(vcf), path(vcf_index)
    path filter_list

    output:
    tuple path("${vcf.getBaseName(2)}_ingroup.vcf.gz"), path("${vcf.getBaseName(2)}_ingroup.vcf.gz.csi")

    script:
    """
    bcftools view --threads 4 -S ^"${filter_list}" --force-samples "${vcf}" -Oz -o "${vcf.getBaseName(2)}_ingroup.vcf.gz"
    bcftools index --threads ${task.cpus} "${vcf.getBaseName(2)}_ingroup.vcf.gz"
    """
}

process pca {
    cpus 4

    input:
    tuple path(vcf), path(vcf_index)

    output:
    path "${vcf.getBaseName(2)}.pca.evec"

    script:
    """
    run_pca.py \
      --vcf "${vcf}" \
      --outprefix "${vcf.getBaseName(2)}" \
      --cleanup \
      --threads ${task.cpus}
    """
}

process bionj {
    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_bionj.nwk"), path("${vcf.getBaseName(2)}_bionj.pdf")

    script:
    """
    run_bionj.py \
      --vcf "${vcf}" \
      --outprefix "${vcf.getBaseName(2)}" \
      --cleanup \
      --r-script run_bionj.R
    """
}

workflow {
    /** Set up channels for input files */
    input_vcf_ch = Channel.fromPath(params.vcfFile, checkIfExists: true)
    bad_samples = Channel.fromPath(params.sampleList.split(',') as List<String>, checkIfExists: true)
    outgroup_samples = Channel.fromPath(params.outgroupList.split(',') as List<String>, checkIfExists: true)
    
    /** Concatenate lists of samples to be filtered from VCFs */
    filter_list = concat_filter_lists(bad_samples.collect())
    outgroup_list = concat_outgroup_lists(outgroup_samples.collect())

    /** Set up filtered VCF files for analysis:
        - a transversions-only VCF -> used for f-statistics
        - a linkage-pruned transversions-only VCF -> used for PCA and ADMIXTURE
        - the pruned, tv-only VCF, restricted to ingroup samples -> used for PCA
    */
    vcf = remove_bad_samples(input_vcf_ch, filter_list)
    transversions = select_transversion_sites(vcf)
    
    pruned = linkage_pruning(transversions)
    ingroup = remove_outgroups(pruned, outgroup_list)

    /** Run analyses */
    pca(ingroup)
    bionj(pruned)

    publish:
    pruned.out to: "${params.outDir}/pruned"
    transversions.out to: "${params.outDir}/transversions"
    pca.out to: "${params.outDir}/pca"
    bionj.out to: "${params.outDir}/bionj"
}
