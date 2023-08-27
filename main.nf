nextflow.enable.dsl=2

/** Input VCF file containing all data */
params.vcfFile

/** Comma-separated list of files listing samples to remove from the VCF */
params.sampleList

/** Files naming outgroup samples to remove from PCA */
params.outgroupList

/** Single file listing the samples and populations to be used for ADMIXTURE */
params.admixtureList

/** Single file listing the population each sample belongs to */
params.populations

/** Single file stating the samples to include in a PCA plot */
params.pcaIngroupSamples

/** A pair of files identifying the samples to exclude from the f4 statistics plots */
params.f4plotExclusions

/** A file that describes how samples are pooled for qpAdm */
params.qpAdmPooling

/** A file listing all the targets to use when computing qpAdm */
params.qpAdmTargets

/** A file listing the population groupings to use when plotting pooled f4 results */
params.f4plotGroups

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
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

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

process annotate_vcf {
    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    input:
    path vcf

    output:
    path "${vcf.getBaseName(2)}_annotated.vcf.gz"

    script:
    """
    make_annotation_bed.sh "${vcf}" # creates id_annotation.bed.gz
    bcftools annotate --threads ${task.cpus} -c CHROM,FROM,TO,ID -a id_annotation.bed.gz -o "${vcf.getBaseName(2)}_annotated.vcf.gz" "${vcf}"
    """
}

process select_transversion_sites {
    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

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
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    publishDir "${params.outDir}/pruned"

    input:
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_pruned.vcf.gz"), path("${vcf.getBaseName(2)}_pruned.vcf.gz.csi")

    script:
    """
    plink --threads ${task.cpus} \
      --indep-pairwise 100kb 1 0.8 \
      --const-fid \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf "${vcf}" \
      --out plink_pruning
    plink --extract plink_pruning.prune.in \
      --const-fid \
      --chr-set 38 \
      --real-ref-alleles \
      --keep-allele-order \
      --vcf "${vcf}" \
      --recode vcf \
      --out "${vcf.getBaseName(2)}_pruned"
    bcftools query -l "${vcf.getBaseName(2)}_pruned.vcf" | sed 's/0_//' > reheader.txt &&
      bcftools reheader --threads ${task.cpus} -s reheader.txt -o tmp "${vcf.getBaseName(2)}_pruned.vcf" &&
      bcftools view --threads ${task.cpus} -o "${vcf.getBaseName(2)}_pruned.vcf.gz" -Oz tmp && rm tmp
    bcftools index --threads ${task.cpus} "${vcf.getBaseName(2)}_pruned.vcf.gz"
    """
}

process remove_outgroups {
    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

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
    // cpus 4
    // executor 'lsf'
    // queue 'normal'
    // memory 2.GB
    // time 2.h
    // clusterOptions '-R "select[avx2]"'

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

process extract_variant_ids {
    input:
    tuple path(vcf), path(vcf_index)

    output:
    path("${vcf.getBaseName(2)}_variant_ids.txt")

    script:
    """
    extract_vcf_variant_ids.py -i "${vcf}" -o "${vcf.getBaseName(2)}_variant_ids.txt"
    """
}

process make_bootstrap_vcf {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 500.MB
    time { 1.h * 2**(task.attempt-1) }
    errorStrategy 'retry'
    maxRetries 2

    input:
    tuple path(vcf), path(vcf_index), path(variant_ids), val(bootstrap_id)

    output:
    tuple path("${vcf.getBaseName(2)}_bootstrap_${bootstrap_id}.vcf.gz"), path("${vcf.getBaseName(2)}_bootstrap_${bootstrap_id}.vcf.gz.csi")

    script:
    """
    make_bootstrap_vcf.py \
      --input-ids="${variant_ids}" \
      --input-vcf="${vcf}" \
      --output-vcf="${vcf.getBaseName(2)}_bootstrap_${bootstrap_id}.vcf.gz" \
      --seed "${bootstrap_id}"
    """
}

process bionj_bootstrap {
    cpus 1
    executor 'lsf'
    queue 'normal'
    memory 500.MB
    time { 1.h * 2**(task.attempt-1) }
    errorStrategy 'retry'
    maxRetries 2

    publishDir "${params.outDir}/bionj/bootstrap"

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

process annotate_bootstrap_info {
    publishDir "${params.outDir}/bionj"

    input:
    path(bionj)
    path('bootstrap/*.nwk')

    output:
    path("${bionj.baseName}.annotated.nwk")

    script:
    """
    bootstrap_annotate_tree.R "${bionj}" bootstrap
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
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

    input:
    path sample_list
    tuple path(vcf), path(vcf_index)

    output:
    tuple path("${vcf.getBaseName(2)}_admixture.bed"), path("${vcf.getBaseName(2)}_admixture.bim"), path("${vcf.getBaseName(2)}_admixture.fam")

    script:
    """
    cut -f1 "${sample_list}" > sample_list.txt
    bcftools view --threads ${task.cpus} -S sample_list.txt --force-samples "${vcf}" -Oz -o "${vcf.getBaseName(2)}_admixture.vcf.gz"
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

process get_admixture_errors {
    input:
    tuple val(k), val(rep), path(logFile)

    output:
    path("${k}_${rep}_error.txt")

    script:
    """
    perl -ne 'if (/^CV error \\(K=(\\d+)\\):\\s+([\\d.]+)\$/) { print "$rep\\t\$1\\t\$2\\n"; }' "$logFile" > err
    perl -ne 'if (/^Loglikelihood: ([\\d.+-]+)\$/) { print "\$1\\n"; }' "$logFile" > ll
    paste err ll > "${k}_${rep}_error.txt" && rm err ll
    """
}

process collate_admixture_errors {
    input:
    path errorFiles

    output:
    path("admixture_errors.txt")

    publishDir "${params.outDir}/admixture"

    script:
    """
    printf "rep\\tk\\terror\\tloglik\\n" > "admixture_errors.txt"
    cat ${errorFiles} | sort -V >> "admixture_errors.txt"
    """
}

process plot_admixture_errors {
    input:
    path errorFile

    output:
    path("${errorFile.baseName}.pdf")

    publishDir "${params.outDir}/admixture"

    script:
    """
    plot_admixture_errors.R \
      -i "${errorFile}" \
      -o "${errorFile.baseName}.pdf"
    """
}

process plot_admixture_result {
    input:
    tuple val(k), val(rep), path(Q), path(P), path(populations)

    output:
    path("${Q.baseName}.pdf")

    publishDir "${params.outDir}/admixture/${rep}"

    script:
    """
    plot_admixture_result.R \
      --input_file "${Q}" \
      --pops "${populations}" \
      --output_file "${Q.baseName}.pdf"
    """
}

process prepare_admixtools_input {
    cpus 4
    executor 'lsf'
    queue 'normal'
    memory 1.GB
    time 2.h

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
    path populations

    output:
    tuple path("${bed.baseName}_HT_combined.pdf"), path("${bed.baseName}_CTVT_combined.pdf")

    publishDir "${params.outDir}/f4stats"

    script:
    """
    f4_statistics.R \
      --input "${bed.baseName}" \
      --output "${bed.baseName}" \
      --populations-file "${populations}" \
      --plot-width 18 \
      --filter "${filter_list[0]}" \
      --second-filter "${filter_list[1]}"
    """
}

process make_qpadm_input_files {
    input:
    tuple path(bed), path(bim), path(fam)
    path(populations)

    output:
    tuple path("${bed.baseName}.pop1.bed"), path("${bed.baseName}.pop1.bim"), path("${bed.baseName}.pop1.fam")
    path("adm_combinations_core.txt")
    path("adm_combinations_comprehensive.txt")

    script:
    """
    repopulate_plink.R \
      --input-prefix "${bed.baseName}" \
      --output-prefix "${bed.baseName}.pop1" \
      --pooling "${populations}"
    make_adm_combinations.R
    """
}

process make_pooled_f4_input_files {
    input:
    tuple path(bed), path(bim), path(fam)
    path(populations)

    output:
    tuple path("${bed.baseName}.pooledf4.bed"), path("${bed.baseName}.pooledf4.bim"), path("${bed.baseName}.pooledf4.fam")

    script:
    """
    repopulate_plink.R \
      --input-prefix "${bed.baseName}" \
      --output-prefix "${bed.baseName}.pooledf4" \
      --pooling "${populations}"
    """
}

process pooled_f4_stats {
    input:
    tuple path(bed), path(bim), path(fam), val(subset)

    output:
    path("${bed.baseName}.${subset}.csv")

    publishDir "${params.outDir}/pooled_f4stats"

    script:
    """
    pooled_f4_statistics.R \
      --input "${bed.baseName}" \
      --output "${bed.baseName}"."${subset}".csv \
      --subset "${subset}"
    """
}

process determine_plotting_order {
    input:
    tuple path(f4result), path(groupings)

    output:
    path("pooled_f4_stats_plotting_order.txt")

    script:
    """
    determine_plotting_order.R \
      --input "${f4result}" \
      --output "pooled_f4_stats_plotting_order.txt" \
      --groupings "${groupings}"
    """
}

process plot_pooled_f4_stats {
    input:
    tuple path(f4result), path(groupings), path(order)

    output:
    tuple path("${f4result.baseName}_ht.pdf"), path("${f4result.baseName}_ctvt.pdf")

    publishDir "${params.outDir}/pooled_f4stats"

    script:
    """
    plot_pooled_f4_statistics.R \
      --input "${f4result}" \
      --output "${f4result}" \
      --order "${order}" \
      --groupings "${groupings}"
    """
}

process run_qpadm {
    cpus { 1 }
    executor 'lsf'
    queue 'small'
    memory 500.MB
    time { 30.min }

    input:
    tuple path(bed), path(bim), path(fam), val(left), val(right), val(target)

    output:
    path("${target}.${left.replace(',','_')}.${right.replace(',','_')}.RDS")

    publishDir "${params.outDir}/qpadm/rds_files"

    script:
    """
    run_adm_combo.R \
      --target="${target}" \
      --left="${left}" \
      --right="${right}" \
      --sampleset="${bed.baseName}"
    """
}

process collate_qpadm {
    input:
    path rds

    output:
    path("qpadm_results.tsv")

    publishDir "${params.outDir}/qpadm"

    script:
    """
    collate_qpadm_results.R \
      --subdir "\$PWD" \
      --outfile "qpadm_results.tsv"
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
    f4plot_exclusions = Channel.fromFilePairs(params.f4plotExclusions, checkIfExists: true)
    qpadm_pooling = Channel.fromPath(params.qpAdmPooling, checkIfExists: true)
    qpadm_targets = Channel.fromPath(params.qpAdmTargets, checkIfExists: true)
    f4plot_groups = Channel.fromPath(params.f4plotGroups, checkIfExists: true)
    
    /** Concatenate lists of samples to be filtered from VCFs */
    filter_list = concat_filter_lists(bad_samples.collect())
    outgroup_list = concat_outgroup_lists(outgroup_samples.collect())

    /** Set up filtered VCF files for analysis:
        - a transversions-only VCF -> used for f-statistics
        - a linkage-pruned transversions-only VCF -> used for PCA and ADMIXTURE
        - the pruned, tv-only VCF, restricted to ingroup samples -> used for PCA
    */
    vcf = remove_bad_samples(input_vcf_ch, filter_list)
    annotated_vcf = annotate_vcf(vcf)
    transversions = select_transversion_sites(annotated_vcf)
    
    pruned = linkage_pruning(transversions)
    ingroup = remove_outgroups(pruned, outgroup_list)

    /** Run analyses */
    /** PCA */
    pca(ingroup)
    plot_pca(pca.out[0], populations, pca_plot_samples)

    /** Phylogeny */
    bionj(pruned)
    variant_ids = extract_variant_ids(pruned)
    /** Run bootstrap replicates in batches of 10 using buffer */
    pruned.combine(variant_ids).combine(Channel.from(1..100)) 
        | make_bootstrap_vcf
        | bionj_bootstrap
    bootstrapped_tree = annotate_bootstrap_info(bionj.out, bionj_bootstrap.out.collect())
    plot_bionj_tree(bootstrapped_tree, populations)
    
    /** Run ADMIXTURE */
    admixture_files = prepare_admixture_input(admixture_samples, pruned)
    admixture_ks = Channel.from(1..15)
    admixture_reps = Channel.from(1..10)
    admixture_inputs = admixture_files.combine(admixture_ks.combine(admixture_reps))
    run_admixture(admixture_inputs)
    run_admixture.out.logs |
      get_admixture_errors |
      collect |
      collate_admixture_errors |
      plot_admixture_errors
    run_admixture.out.results.combine(admixture_samples) |
      plot_admixture_result

    /** f4-statistics */
    admixtools_files = prepare_admixtools_input(transversions)
    f4_stats(admixtools_files, f4plot_exclusions, populations)

    /** pooled f4-statistics */
    make_pooled_f4_input_files(admixtools_files, populations)
        | combine(Channel.of("all", "chr1", "chr7", "chr21"))
        | pooled_f4_stats

    /** f4-statistics plots */
    plotting_order_input = pooled_f4_stats.out |
      branch { it ->
          TRUE: it.baseName.endsWith("all")
      }
    plotting_order = determine_plotting_order(plotting_order_input.TRUE.combine(f4plot_groups))
    plot_pooled_f4_stats(pooled_f4_stats.out.combine(f4plot_groups).combine(plotting_order)) 

    /** qpAdm */
    qpadm_inputs = make_qpadm_input_files(admixtools_files, qpadm_pooling)
    combos = qpadm_inputs[1] | splitCsv(sep: '\t', header: true)
    targets = qpadm_targets | splitText() { it.trim() }
    run_qpadm_input_ch = qpadm_inputs[0].combine(combos).combine(targets)
      .map { it -> [it[0], it[1], it[2], it[3].LEFT, it[3].RIGHT, it[4] ]}
      .filter { 
        def left = it[3].split(',');
        def right = it[4].split(',');
        def target = it[5];
        !(left.contains(target) || right.contains(target))
      }
    qpadm_results = run_qpadm(run_qpadm_input_ch)
    qpadm_results.collect() | collate_qpadm
}
