nextflow.enable.dsl=2

/** Input VCF file containing all data */
params.vcfFile

/** Comma-separated list of files listing samples to remove from the VCF */
params.sampleList

/** Files naming outgroup samples to remove from PCA */
params.outgroupList

/** Single file listing the samples to be used for ADMIXTURE */
params.admixtureList

/** Single file listing the population each sample belongs to */
params.populations

/** Single file stating the samples to include in a PCA plot */
params.pcaIngroupSamples

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
    publishDir "${params.outDir}/transversions"

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
    publishDir "${params.outDir}/pruned"

    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_pruned.vcf.gz"), path("${vcf.getBaseName(2)}_pruned.vcf.gz.csi")

    script:
    """
    make_annotation_bed.sh "${vcf}" # creates id_annotation.bed.gz
    bcftools annotate --threads ${task.cpus} -c CHROM,FROM,TO,ID -a id_annotation.bed.gz -o annotated.vcf.gz "${vcf}"
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
    rm annotated.vcf.gz
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
    publishDir "${params.outDir}/pca"

    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}.pca.evec"), path("${vcf.getBaseName(2)}.eval")

    script:
    """
    run_pca.py \
      --vcf "${vcf}" \
      --outprefix "${vcf.getBaseName(2)}" \
      --cleanup \
      --threads ${task.cpus}
    """
}

process plot_pca {
    publishDir "${params.outDir}/pca"

    input:
    tuple path(pca_evec_file), path(pca_eval_file)
    path populations
    path samples

    output:
    path "${pca_evec_file.baseName}.pdf"

    script:
    """
    pca_plot.R -p "${populations}" -i "${samples}" -d "${pca_evec_file}" -o "${pca_evec_file.baseName}.pdf"
    """
}

process bionj {
    publishDir "${params.outDir}/bionj"

    input:
    tuple path(vcf), path(vcf_index)

    output:
    path("${vcf.getBaseName(2)}_bionj.nwk")

    script:
    """
    run_bionj.py \
      --vcf "${vcf}" \
      --outprefix "${vcf.getBaseName(2)}" \
      --cleanup \
      --r-script run_bionj.R
    """
}

process plot_bionj_tree {
    publishDir "${params.outDir}/bionj"

    input:
    path(newick)
    path(populations)

    output:
    path("${newick.baseName}.pdf")

    script:
    """
    tree_plot.R -p "${populations}" -t "${newick}" -o "${newick.baseName}.pdf"
    """
}

process prepare_admixture_input {
    cpus 4

    input:
    path sample_list
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_admixture.bed"), path("${vcf.getBaseName(2)}_admixture.bim"), path("${vcf.getBaseName(2)}_admixture.fam")

    script:
    """
    bcftools view --threads ${task.cpus} -S "${sample_list}" --force-samples "${vcf}" -Oz -o "${vcf.getBaseName(2)}_admixture.vcf.gz"
    plink --threads ${task.cpus} \
      --double-id \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf "${vcf.getBaseName(2)}_admixture.vcf.gz" \
      --make-bed \
      --out "${vcf.getBaseName(2)}_admixture"
    """
}

process run_admixture {
    cpus { 4 * 2**(task.attempt-1) }
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time { 2.h * 2**(task.attempt-1) }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.outDir}/admixture/${rep}"

    input:
    tuple path(bed), path(bim), path(fam), val(k), val(rep)

    output:
    tuple val(k), val(rep), path("${bed.baseName}.${k}.Q"), path("${bed.baseName}.${k}.P"), emit: results
    tuple val(k), val(rep), path("${bed.baseName}.${k}.log"), emit: logs

    script:
    seed = "99${k * 1000}${rep}".toString()
    
    """
    admixture --cv=10 -j${task.cpus} -s "${seed}" "${bed}" ${k} | tee "${bed.baseName}.${k}.log"
    """
}

process prepare_admixtools_input {
    cpus 4

    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_admixtools.bed"), path("${vcf.getBaseName(2)}_admixtools.bim"), path("${vcf.getBaseName(2)}_admixtools.fam")

    script:
    """
    plink --threads ${task.cpus} \
      --double-id \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf "${vcf}" \
      --make-bed \
      --out "${vcf.getBaseName(2)}_admixtools"
    """
}
process f4_stats {
    input:
    tuple path(bed), path(bim), path(fam)
    tuple val(id), file(filter_list)

    output:
    path("${vcf.getBaseName(2)}_f4_stats.pdf")

    script:
    """
    f4_statistics.R \
      --input "${bed}" \
      --output "${vcf.baseName}" \
      --plot-width 18 \
      --filter "${filter_list[0]}" \
      --second-filter "${filter_list[1]}"
    """
}

workflow {
    /** Set up channels for input files */
    input_vcf_ch = Channel.fromPath(params.vcfFile, checkIfExists: true)
    bad_samples = Channel.fromPath(params.sampleList.split(',') as List<String>, checkIfExists: true)
    outgroup_samples = Channel.fromPath(params.outgroupList.split(',') as List<String>, checkIfExists: true)
    admixture_samples = Channel.fromPath(params.admixtureList, checkIfExists: true)
    populations = Channel.fromPath(params.populations, checkIfExists: true)
    pca_plot_samples = Channel.fromPath(params.pcaIngroupSamples, checkIfExists: true)
    f4plot_exclusions = Channel.fromFilePairs("data/sample_lists/f4_exclusions_{1,2}.txt", checkIfExists: true)
    
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
    /** PCA */
    pca(ingroup)
    plot_pca(pca.out[0], populations, pca_plot_samples)

    /** Phylogeny */
    bionj(pruned)
    plot_bionj_tree(bionj.out[0], populations)
    
    /** Run ADMIXTURE */
    admixture_files = prepare_admixture_input(admixture_samples, pruned)
    admixture_ks = Channel.from(1..15)
    admixture_reps = Channel.from(1..10)
    admixture_inputs = admixture_files.combine(admixture_ks.combine(admixture_reps))
    run_admixture(admixture_inputs)

    /** f4-statistics */

    admixtools_files = prepare_admixtools_input(transversions)
    f4_stats(transversions, f4plot_exclusions)
}
