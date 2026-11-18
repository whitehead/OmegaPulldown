/*
 * Pair construction (predictor-agnostic) and Boltz-2 YAML assembly.
 *
 * GENERATE_PAIRS resolves the run's two inputs, one pool
 * (params.sequences_fasta) and one pair table (params.pairs_tsv), into the
 * pairs.tsv that drives every downstream stage.
 *
 * This process is the resolver between a permissive input schema (first two
 * columns are IDs, positional, header auto-detected: so a raw proteomics
 * export drops in unmodified) and a strict emitted one (pair_id, bait_id,
 * bait_seq, candidate_id, candidate_seq). Those two schemas are deliberately
 * NOT the same thing: the input is an edge list, the emitted file is the
 * materialized screen with sequences inlined.
 *
 * CHAIN ORDER: column 1 is the bait is chain A, guaranteed from here on.
 * Pairs are collapsed on tuple(sorted((a,b))) since A_with_B and B_with_A are
 * the same complex, and the FIRST orientation listed is the one kept: so
 * the collapse discards a later restatement rather than re-orienting the row
 * that survives.
 *
 * pairs.tsv is also where the MSA universe comes from (main.nf derives the
 * unique-protein channel from this file's rows, not from any FASTA).
 *
 * MAKE_BOLTZ_YAML / MAKE_BOLTZ_YAML_NO_MSA then assemble one Boltz-2 input
 * YAML per pair, with or without pre-computed MSA CSVs. (The AF3 equivalent,
 * MAKE_AF3_JSON, lives in predict_af3.nf.)
 *
 * PTMs come from ONE id-keyed table (params.ptms, conventionally ptms.tsv;
 * tab or comma). make_boltz_yaml.py resolves each chain's modifications by
 * protein ID, never by whether that protein is the bait or the candidate of
 * this particular pair. See resolvePtms() in main.nf.
 */

process GENERATE_PAIRS {
    label 'process_low'

    publishDir "${params.outdir}/pairs", mode: 'copy'

    input:
    path pairs_tsv,       stageAs: 'pairgen_pairs.tsv'
    path sequences_fasta, stageAs: 'pairgen_sequences.fasta'

    output:
    path "pairs.tsv", emit: pairs_tsv

    script:
    """
    generate_pairs.py \
        --pairs-tsv ${pairs_tsv} \
        --sequences ${sequences_fasta} \
        --output pairs.tsv
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real; also
    // gives downstream stub processes a genuine pairs.tsv to fan out over.
    """
    generate_pairs.py \
        --pairs-tsv ${pairs_tsv} \
        --sequences ${sequences_fasta} \
        --output pairs.tsv
    """
}


process MAKE_BOLTZ_YAML {
    tag "${pair_id}"
    label 'process_low'

    // Content-addressed caching: the Boltz-2 half of the same fix applied
    // to MAKE_AF3_JSON.
    //
    cache 'deep'

    publishDir "${params.outdir}/yamls", mode: 'copy'

    input:
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq), path(bait_msa), path(candidate_msa)
    val ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), path(bait_msa), path(candidate_msa), emit: yaml_with_msa

    script:
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real.
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --bait-msa "${bait_msa}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        --candidate-msa "${candidate_msa}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """
}

/*
 * Alternative: generate YAML without pre-computed MSA.
 * Boltz-2 will fetch MSAs from the server at prediction time.
 * Simpler but recomputes MSAs for shared proteins.
 */
process MAKE_BOLTZ_YAML_NO_MSA {
    tag "${pair_id}"
    label 'process_low'

    publishDir "${params.outdir}/yamls", mode: 'copy'

    input:
    tuple val(pair_id), val(bait_id), val(bait_seq), val(candidate_id), val(candidate_seq)
    val ptms_path

    output:
    tuple val(pair_id), path("${pair_id}.yaml"), emit: yaml

    script:
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """

    stub:
    // Pure-Python, local files only: safe and fast to run for real.
    // One ID-keyed PTM table (or a comma-separated list of them); the
    // script resolves each chain's PTMs by protein ID, never by bait or
    // candidate role.
    def ptms_arg = (ptms_path != 'no_file') ? "--ptms ${ptms_path}" : ""
    """
    make_boltz_yaml.py \
        --bait-id "${bait_id}" \
        --bait-seq "${bait_seq}" \
        --candidate-id "${candidate_id}" \
        --candidate-seq "${candidate_seq}" \
        ${ptms_arg} \
        --output "${pair_id}.yaml"
    """
}
