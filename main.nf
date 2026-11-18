/*
 * =========================================================================
 *  OmegaPulldown: Protein-protein interaction screen using Boltz-2
 *                 or AlphaFold 3
 * =========================================================================
 *
 *  Input:  one protein pool and a table for paired proteins (complexes).
 *            --sequences_fasta   a FASTA holding every protein in the study
 *            --pairs_tsv         the pairs to predict; the first two columns
 *                                are protein IDs, header auto-detected
 *          Initialize with bin/make_pairs.py, which writes an ordinary pair
 *          list.
 *  Output: ranked table of predicted complex structures with confidence
 *          metrics (ipTM, pDockQ, pLDDT, pTM), plus per-pair structure
 *          subdirectories under results/structures_boltz2/ or
 *          results/structures_af3/ depending on the predictor/folding backend.
 *
 *  Pipeline stages:
 *    1. GENERATE_PAIRS - initialize the pair list: resolve each ID in
 *                        --pairs_tsv against --sequences_fasta and inline
 *                        both sequences. Pairs are collapsed on
 *                        tuple(sorted((a,b))), since A_with_B and B_with_A
 *                        are the same complex, keeping the first orientation
 *                        listed. The emitted pairs.tsv is also where the MSA
 *                        universe comes from.
 *
 *                        CHAIN ORDER: column 1 is the bait is chain A, from
 *                        the input list through to the output structures.
 *
 *    2a. (Boltz-2 local MSA, use_msa_server=false)
 *        MSA_SEARCH_DB_BOLTZ2, UniRef30 search per unique protein (CPU);
 *        then (only if params.use_env=true) MSA_SEARCH_ENV_BOLTZ2,
 *        env/metagenomic search; then MERGE_TAXONOMY_BOLTZ2, merge +
 *        taxonomy-ID tagging into the final per-protein CSVs (see
 *        modules/search_msa_boltz2.nf). No separate pairing/join step:
 *        Boltz-2's own featurizer pairs rows sharing a taxonomy ID across
 *        whichever two chains' CSVs are loaded for a prediction.
 *
 *    2b. (AF3-native MSA, predictor=alphafold3)
 *        MSA_SEARCH_DB, one task per unpaired database (UniRef90, MGnify,
 *        Small BFD); then MERGE_AF3_NATIVE_MSA, merge with exact-duplicate
 *        removal into one unpaired a3m per protein. Separately,
 *        MSA_SEARCH_UNIPROT_PAIRED builds the PAIRED a3m from UniProt, which
 *        AF3 keeps distinct from the unpaired MSA and uses for cross-chain
 *        species matching (see modules/search_msa_af3.nf).
 *
 *    3. MAKE_BOLTZ_YAML / MAKE_AF3_JSON - assemble input per pair
 *    4. BOLTZ_PREDICT / AF3_PREDICT     - structure prediction (GPU)
 *    5. PARSE_METRICS / PARSE_METRICS_AF3 - confidence + pDockQ per pair
 *                        (pDockQ via bin/compute_pdockq.py, shared by both)
 *    6. RANK_PREDICTIONS - aggregate and rank all pairs
 *
 *  Predictor selection:
 *    --predictor boltz2       (default)
 *    --predictor alphafold3
 *
 *  MSA caching (note that the two predictors cache different intermediates:
 *    predictor=boltz2 caches per-protein CSVs at
 *      params.alignments_dir_boltz2/csv/<prot_id>.csv
 *    There is no a3m-level cache for Boltz-2: taxonomy tagging happens in
 *    the same task as the search, so no intermediate a3m is worth keeping.
 *    predictor=alphafold3 caches per-protein a3m files at
 *      params.alignments_dir_af3/a3m/<prot_id>.a3m
 *      params.alignments_dir_af3/paired_a3m/<prot_id>.paired.a3m
 *    A protein whose footprint already exists skips its search entirely.
 *    Boltz-2 and AF3 use separate alignment directories: they were trained
 *    on different MSA lineages, so a shared directory is not a neutral cache.
 *
 *  When params.use_msa_server = true, stage 2a is skipped and Boltz-2
 *  fetches MSAs from the ColabFold API at prediction time.
 *  (Only valid with predictor = boltz2.)
 *
 *  Post-translational modifications (optional):
 *    Place ptms.tsv in the same directory as --sequences_fasta, or specify
 *    it explicitly via --ptms. Modifications are resolved by protein ID, so
 *    one table covers both chains of every pair.  Format (tab or comma):
 *      protein_id	position	ccd
 *      BAIT1_GFP	18	SEP
 *
 *  Usage:
 *    All runs take the same two inputs. If you want a cross product or an
 *    all-vs-all rather than an explicit list, generate the list first:
 *
 *      bin/make_pairs.py --all-vs-all pool.fasta > pairs/all_vs_all.tsv
 *      bin/make_pairs.py --baits baits.fasta --candidates candidates.fasta \
 *                        --sequences pool.fasta > pairs/pulldown.tsv
 *
 *    The extra step is deliberate: the pair set becomes a table that one can
 *    read, count and diff before committing GPU time to it.
 *
 *    # Boltz-2 (default), API-based MSA (use_msa_server=true is the default:
 *    # recommended unless running a large all-vs-all screen, see
 *    # boltz2_uniref30_db/taxonomy_rdb in nextflow.config for the local-MSA
 *    # alternative)
 *    nextflow run main.nf \
 *      --sequences_fasta inputs/sequences/proteins.fasta \
 *      --pairs_tsv inputs/pairs/actin_pulldown.tsv
 *
 *    # Boltz-2, local MSA mode (no public API dependency: recommended for
 *    # large all-vs-all screens; see modules/search_msa_boltz2.nf). UniRef30
 *    # search only by default (params.use_env=false, see nextflow.config
 *    # for why); add --use_env true --boltz2_envdb_db /path/to/envdb to
 *    # also run the env/metagenomic search, at the cost of a scarce
 *    # bigmem-class node allocation.
 *    nextflow run main.nf \
 *      --sequences_fasta inputs/sequences/proteins.fasta \
 *      --pairs_tsv inputs/pairs/all_vs_all.tsv \
 *      --use_msa_server false \
 *      --boltz2_uniref30_db /path/to/uniref30_2302_db \
 *      --taxonomy_rdb /path/to/taxonomy.rdb
 *
 *    # AlphaFold 3 (use_msa_server is not supported here: AF3 always uses
 *    # the locally built AF3-native MSAs of stage 2b)
 *    nextflow run main.nf \
 *      --predictor alphafold3 \
 *      --use_msa_server false \
 *      --sequences_fasta inputs/sequences/proteins.fasta \
 *      --pairs_tsv inputs/pairs/actin_pulldown.tsv
 *
 *    # An experimental edge list straight from proteomics: predict ONLY the
 *    # pairs reported, to add a structure-based ranking to evidence that
 *    # already exists. No rule to evaluate: the list IS the pair set. The
 *    # first two columns are protein IDs (header auto-detected) and further
 *    # columns are ignored, so a SAINT/IP-MS export usually drops in
 *    # unmodified.
 *    nextflow run main.nf \
 *      --sequences_fasta all_proteins.fasta \
 *      --pairs_tsv ms_hits.tsv
 *
 *    # MSA stages only, no GPU work (--msa_only). For validating MSA
 *    # infrastructure (a moved reference DB, a rebuilt index, a new node
 *    # type) without competing for GPU partitions that real screens are
 *    # using. The per-protein cache is still written, so a later run without
 *    # the flag reuses it and goes straight to prediction.
 *    nextflow run main.nf \
 *      --sequences_fasta all_proteins.fasta \
 *      --pairs_tsv ms_hits.tsv \
 *      --use_msa_server false \
 *      --msa_only true
 *
 *    # Reuse pre-computed alignments (skips the MSA stage for any protein
 *    # whose cached artefact is already present, see "MSA caching" above)
 *    nextflow run main.nf \
 *      --predictor alphafold3 \
 *      --use_msa_server false \
 *      --sequences_fasta all_proteins.fasta \
 *      --pairs_tsv ms_hits.tsv \
 *      --alignments_dir_af3 results/alignments_af3
 * =========================================================================
 */

nextflow.enable.dsl = 2

// --- Module imports ---
include { MSA_SEARCH_DB_BOLTZ2; MSA_SEARCH_ENV_BOLTZ2; MERGE_TAXONOMY_BOLTZ2 }                 from './modules/search_msa_boltz2'
include { MSA_SEARCH_DB; MERGE_AF3_NATIVE_MSA; MSA_SEARCH_UNIPROT_PAIRED }                     from './modules/search_msa_af3'
include { GENERATE_PAIRS; MAKE_BOLTZ_YAML; MAKE_BOLTZ_YAML_NO_MSA }   from './modules/prepare_input'
include { BOLTZ_PREDICT }                                              from './modules/predict_boltz2'
include { PARSE_METRICS; RANK_PREDICTIONS }                            from './modules/rank'
include { MAKE_AF3_JSON; AF3_PREDICT }                                 from './modules/predict_af3'
include { PARSE_METRICS_AF3 }                                          from './modules/rank_af3'


/*
 * Resolve the PTM table for this run.
 *
 * ONE ID-keyed table, conventionally ptms.tsv:
 *   protein_id<TAB>position<TAB>ccd
 * Tab- or comma-separated (checked). Resolution by protein ID alone, never
 * by whether a protein is the bait or the candidate of a given pair.
 *
 * Priority:
 *   1. Explicit --ptms
 *   2. Convention: ptms.tsv, then ptms.csv, beside sequences_fasta
 *   3. No PTMs ('no_file' sentinel)
 */
def resolvePtms(poolFasta) {
    // 1. Explicit
    if (params.ptms != null) {
        def f = file(params.ptms)
        if (!f.exists()) {
            error "ptms file not found: ${params.ptms}"
        }
        log.info "PTM table (explicit): ${f}"
        return f.toString()
    }

    // 2. Convention: a sibling of the pool. PTMs annotate PROTEINS, not
    // screens (a protein's phospho-sites do not change because it was
    // paired differently) so the table belongs beside the sequences and is
    // reused by every pair list run against that pool.
    if (poolFasta != null) {
        def dir = file(poolFasta).parent
        for (name in ['ptms.tsv', 'ptms.csv']) {
            def p = dir.resolve(name)
            if (p.exists()) {
                log.info "PTM table (auto-detected): ${p}"
                return p.toString()
            }
        }
    }

    // 3. None
    return 'no_file'
}

workflow {

    // ---------------------------------------------------------------
    // Validate predictor choice
    // ---------------------------------------------------------------
    def predictor = params.predictor.toLowerCase()
    if (!(predictor in ['boltz2', 'alphafold3'])) {
        error "Invalid predictor '${params.predictor}'. Choose 'boltz2' or 'alphafold3'."
    }
    if (predictor == 'alphafold3' && params.use_msa_server) {
        error "use_msa_server is not supported with predictor=alphafold3 " +
              "(it fetches Boltz-2-lineage MSAs from the ColabFold API; AF3 " +
              "needs its own AF3-native MSAs). Pass --use_msa_server false. " +
              "Note -profile test sets it true, so that profile needs the " +
              "override explicitly."
    }

    // msa_only: stop once the per-protein MSA output exists, before any
    // structure predictions are submitted.
    if (params.msa_only && params.use_msa_server && predictor == 'boltz2') {
        error "--msa_only has nothing to do with --use_msa_server true: " +
              "in server mode Boltz-2 fetches its MSAs during prediction, so " +
              "there is no MSA stage to stop after. Pass --use_msa_server false."
    }

    // ---------------------------------------------------------------
    // Validate the input model: single protein pool, single pair table
    // ---------------------------------------------------------------

    // List everything that is missing, rather than failing on the first one.
    // All input params default to null.
    def missing = []
    if (!params.sequences_fasta)
        missing << '--sequences_fasta  a FASTA holding every protein the run touches'
    if (!params.pairs_tsv)
        missing << '--pairs_tsv        the pairs to predict (first two columns are protein IDs)'
    if (missing) {
        error "Missing required input(s):\n" +
              missing.collect { "    ${it}" }.join('\n') +
              "\n\nTo build a pair list from a rule rather than writing one:\n" +
              "    bin/make_pairs.py --all-vs-all pool.fasta > pairs.tsv\n" +
              "    bin/make_pairs.py --baits b.fasta --candidates c.fasta \\\n" +
              "                      --sequences pool.fasta > pairs.tsv"
    }

    def pool_file  = file(params.sequences_fasta)
    def pairs_file = file(params.pairs_tsv)
    if (!pool_file.exists())
        error "sequences_fasta not found: ${params.sequences_fasta}"
    if (!pairs_file.exists())
        error "pairs_tsv not found: ${params.pairs_tsv}"

    // Detect a bad --rank_by up front, not at the last process of the run.
    def shared_metrics = ['ptm', 'iptm', 'protein_iptm',
                          'pdockq', 'n_interface_contacts', 'avg_interface_plddt']
    def af3_metrics    = ['af3_ranking_score', 'af3_complex_plddt', 'af3_complex_iplddt']
    def boltz2_metrics = ['boltz2_confidence_score', 'boltz2_complex_plddt',
                          'boltz2_complex_iplddt']
    def is_af3         = params.predictor.toLowerCase() == 'alphafold3'
    def known_metrics  = shared_metrics + (is_af3 ? af3_metrics : boltz2_metrics)
    def other_metrics  = is_af3 ? boltz2_metrics : af3_metrics

    def rb = (params.rank_by ?: '').replaceAll(/_mean$/, '')
    if (!known_metrics.contains(rb)) {
        // Name the specific mistake when it is the other predictor's column.
        //
        // Indent the hint to match the block below: stripIndent() removes
        // the common leading whitespace, so a flush-left interpolated line
        // would silently defeat it for the whole message.
        def hint = other_metrics.contains(rb)
            ? "\n            That is the other predictor's column; this run is " +
              "predictor=${params.predictor}.\n"
            : ''
        error """
            --rank_by '${params.rank_by}' is not a valid metric for predictor=${params.predictor}.
            ${hint}
            Comparable across both predictors:
                ${shared_metrics.join(', ')}

            Specific to this predictor:
                ${(is_af3 ? af3_metrics : boltz2_metrics).join(', ')}

            A bare name is shared and may be read across both rankings
            tables; an af3_ / boltz2_ prefix marks a number that is only
            meaningful within its own table. See the column reference in
            modules/rank_af3.nf.
            """.stripIndent()
    }

    // Backend-specific alignments directory: Boltz-2 and AF3 were trained
    // on different MSA lineages, so they don't share a cache.
    def alignments_dir = (predictor == 'alphafold3') ? params.alignments_dir_af3 : params.alignments_dir_boltz2

    log.info """
    =========================================================================
    OmegaPulldown
    =========================================================================
    Predictor    : ${predictor}
    Pool         : ${pool_file}
    Pairs        : ${pairs_file}
    Alignments   : ${alignments_dir}
    MSA server   : ${params.use_msa_server}
    Output       : ${params.outdir}
    =========================================================================
    """.stripIndent()

    // ---------------------------------------------------------------
    // Resolve optional PTM CSV files
    // ---------------------------------------------------------------
    // One pool, so no nomination to make: the PTM table is --ptms, or a
    // sibling of sequences_fasta, or nothing.
    ptms_path = resolvePtms(params.sequences_fasta)

    // ---------------------------------------------------------------
    // Stage 1: Create the pair list
    //
    // One pool (sequences_fasta) and one pair table (pairs_tsv), resolved
    // into pairs.tsv with both sequences inlined.
    // ---------------------------------------------------------------
    pairs_input_ch = Channel.fromPath(params.pairs_tsv)
    pool_ch        = Channel.fromPath(params.sequences_fasta)

    GENERATE_PAIRS(pairs_input_ch, pool_ch)

    // Parse TSV -> channel of (pair_id, bait_id, bait_seq, candidate_id, candidate_seq)
    pairs_ch = GENERATE_PAIRS.out.pairs_tsv
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(
            row.pair_id,
            row.bait_id,
            row.bait_seq,
            row.candidate_id,
            row.candidate_seq
        )}

    if (params.use_msa_server && predictor == 'boltz2') {

        // =============================================================
        // Server mode: skip MSA pre-computation (Boltz-2 only)
        // Boltz-2 will use --use_msa_server to fetch MSAs on the fly
        // =============================================================

        predict_ch = MAKE_BOLTZ_YAML_NO_MSA(
            pairs_ch,
            ptms_path
        ).yaml
            .map { pair_id, yaml_file -> tuple(pair_id, yaml_file, []) }

    } else {

        // =============================================================
        // Local MSA mode: pre-compute MSAs per unique protein
        // Avoids redundant searches when proteins appear in many pairs.
        // Skips proteins that already have A3M files in alignments_dir.
        // =============================================================

        // The MSA universe is read from pairs.tsv, not from the pool.
        //
        unique_proteins_ch = pairs_ch
            .flatMap { pair_id, bait_id, bait_seq, cand_id, cand_seq ->
                [tuple(bait_id, bait_seq), tuple(cand_id, cand_seq)] }
            .unique { it[0] }

        // Stage 2: MSA search for proteins without cached MSAs (CPU-bound).
        // predictor=alphafold3 uses the AF3-native 4-DB search,
        // caching per-protein a3m files; predictor=boltz2 uses
        // the ColabFold-lineage single-DB search + taxonomy tagging.
        if (predictor == 'alphafold3') {
            // --- MSA caching: split into cached vs. need-to-compute ---
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/a3m/${it[0]}.a3m").exists()
                compute: true
            }.set { split_ch }

            // Use cached A3M files directly
            cached_a3m = split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/a3m/${prot_id}.a3m"))
            }

            // Batch proteins needing search into fixed-size groups before
            // invoking MSA_SEARCH_DB. This happens AFTER the cache-hit
            // branch above, not before, so caching still decides per-protein
            // whether a protein needs work at all.
            db_ch = Channel.from(params.af3_native_dbs)
            // .collate() emits each batch as a List; .combine() flattens
            // List-valued items by one level when building its output tuple.
            batched_compute_ch = split_ch.compute.collate(params.af3_native_batch_size).map { batch -> [batch] }
            search_input_ch = batched_compute_ch.combine(db_ch)

            MSA_SEARCH_DB(search_input_ch)

            // MSA_SEARCH_DB now emits one file-list per (batch, db) task,
            // not one (prot_id, path) tuple per protein: flatten back into
            // individual files, then recover prot_id from each filename
            // (format "<prot_id>.<db.name>.a3m", named via mmseqs2's
            // --unpack-name-mode 0, see MSA_SEARCH_DB in
            // modules/search_msa_af3.nf). Requires prot_id to contain no
            // literal '.'.
            per_protein_a3m_ch = MSA_SEARCH_DB.out.a3m_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }
                .groupTuple()

            MERGE_AF3_NATIVE_MSA(per_protein_a3m_ch)
            computed_a3m = MERGE_AF3_NATIVE_MSA.out.a3m

            // Merge cached + newly computed A3M streams
            all_a3m = cached_a3m.mix(computed_a3m)

            // Paired MSA (UniProt only): independent cache check.
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/paired_a3m/${it[0]}.paired.a3m").exists()
                compute: true
            }.set { paired_split_ch }

            cached_paired_msa = paired_split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/paired_a3m/${prot_id}.paired.a3m"))
            }

            // Same batching as the unpaired search above.
            batched_paired_compute_ch = paired_split_ch.compute.collate(params.af3_native_batch_size)
            MSA_SEARCH_UNIPROT_PAIRED(batched_paired_compute_ch)

            computed_paired_msa = MSA_SEARCH_UNIPROT_PAIRED.out.a3m_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }

            all_paired_msa = cached_paired_msa.mix(computed_paired_msa)
        } else {
            // predictor == 'boltz2', local MSA mode (use_msa_server=false)

            // --- CSV caching: split into cached vs. need-to-compute ---
            // (CSV, not a3m: see the comment above this if/else for why
            // this predictor caches at the CSV level directly.)
            unique_proteins_ch.branch {
                cached:  file("${alignments_dir}/csv/${it[0]}.csv").exists()
                compute: true
            }.set { csv_split_ch }

            cached_csv = csv_split_ch.cached.map { prot_id, seq ->
                tuple(prot_id, file("${alignments_dir}/csv/${prot_id}.csv"))
            }

            // Batch proteins needing search+tagging, same rationale as
            // AF3-native's batching above (and MSA_SEARCH_UNIPROT_PAIRED).
            batched_compute_ch = csv_split_ch.compute.collate(params.boltz2_batch_size)
            MSA_SEARCH_DB_BOLTZ2(batched_compute_ch)

            if (params.use_env) {
                MSA_SEARCH_ENV_BOLTZ2(MSA_SEARCH_DB_BOLTZ2.out.uniref_result)
                merge_input_ch = MSA_SEARCH_ENV_BOLTZ2.out.env_result
            } else {
                no_env_file = file("${projectDir}/assets/NO_ENV")
                merge_input_ch = MSA_SEARCH_DB_BOLTZ2.out.uniref_result
                    .map { batch, qdb, prof_res, uniref_a3m -> tuple(batch, qdb, uniref_a3m, no_env_file) }
            }
            MERGE_TAXONOMY_BOLTZ2(merge_input_ch)

            computed_csv = MERGE_TAXONOMY_BOLTZ2.out.csv_files
                .flatten()
                .map { f -> tuple(f.name.tokenize('.')[0], f) }

            all_csv = cached_csv.mix(computed_csv)
        }

        // --- Branch by predictor ---
        //
        // Stage 4 (assembling one Boltz-2 YAML / AF3 JSON per pair) is
        // skipped under --msa_only too, not just the GPU stages below.
        if (params.msa_only) {
            // Nothing to assemble; the MSA artefacts above are the product.
        }
        else if (predictor == 'boltz2') {

            // Stage 4: Join MSA CSVs with bait-candidate pairs
            pairs_with_bait_msa = pairs_ch
                .map { pair_id, bait_id, bait_seq, candidate_id, candidate_seq ->
                    tuple(bait_id, pair_id, bait_seq, candidate_id, candidate_seq)
                }
                .combine(all_csv, by: 0)
                .map { bait_id, pair_id, bait_seq, candidate_id, candidate_seq, bait_csv ->
                    tuple(candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_csv)
                }

            pairs_with_both_msa = pairs_with_bait_msa
                .combine(all_csv, by: 0)
                .map { candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_csv, candidate_csv ->
                    tuple(pair_id, bait_id, bait_seq, candidate_id, candidate_seq, bait_csv, candidate_csv)
                }

            predict_ch = MAKE_BOLTZ_YAML(
                pairs_with_both_msa,
                ptms_path
            ).yaml_with_msa
                .map { pair_id, yaml, bait_csv, cand_csv -> tuple(pair_id, yaml, [bait_csv, cand_csv]) }

        } else {
            // predictor == 'alphafold3'

            // Stage 4: Join unpaired + paired A3M files with bait-candidate pairs
            pairs_with_bait_msa = pairs_ch
                .map { pair_id, bait_id, bait_seq, candidate_id, candidate_seq ->
                    tuple(bait_id, pair_id, bait_seq, candidate_id, candidate_seq)
                }
                .combine(all_a3m, by: 0)
                .combine(all_paired_msa, by: 0)
                .map { bait_id, pair_id, bait_seq, candidate_id, candidate_seq, bait_a3m, bait_paired_msa ->
                    tuple(candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m, bait_paired_msa)
                }

            pairs_with_both_msa = pairs_with_bait_msa
                .combine(all_a3m, by: 0)
                .combine(all_paired_msa, by: 0)
                .map { candidate_id, pair_id, bait_id, bait_seq, candidate_seq, bait_a3m, bait_paired_msa, candidate_a3m, candidate_paired_msa ->
                    tuple(pair_id, bait_id, bait_seq, candidate_id, candidate_seq, bait_a3m, candidate_a3m, bait_paired_msa, candidate_paired_msa)
                }

            predict_ch = MAKE_AF3_JSON(
                pairs_with_both_msa,
                ptms_path
            ).json_with_msa
                .map { pair_id, json_file, a3m_files -> tuple(pair_id, json_file, a3m_files) }
        }
    }


    // ---------------------------------------------------------------
    // Stages 5-7: prediction, metrics, ranking: all skipped under --msa_only
    // ---------------------------------------------------------------
    if (params.msa_only) {
        log.info "--msa_only: stopping after the MSA stage, before any GPU work."
    }
    else {
        // Stage 5: Structure prediction (GPU-bound)
        if (predictor == 'boltz2') {
            BOLTZ_PREDICT(predict_ch)
            pred_results = BOLTZ_PREDICT.out.results
        } else {
            AF3_PREDICT(predict_ch)
            pred_results = AF3_PREDICT.out.results
        }

        // Stage 6: Parse metrics per pair (pDockQ + confidence)
        if (predictor == 'boltz2') {
            PARSE_METRICS(pred_results)
            all_metrics = PARSE_METRICS.out.metrics.collect()
        } else {
            PARSE_METRICS_AF3(pred_results)
            all_metrics = PARSE_METRICS_AF3.out.metrics.collect()
        }

        // Stage 7: Aggregate and rank all predictions
        RANK_PREDICTIONS(all_metrics)
    }
}


workflow.onComplete {

    // -------------------------------------------------------------------
    // Reconcile pairs enumerated against pairs actually ranked.
    //
    // The process-wide errorStrategy in nextflow.config is
    //   { task.attempt <= 2 ? 'retry' : 'ignore' }
    // and that is deliberate: on a several-hundred-candidate screen, one
    // pair that OOMs or hits a bad node must not take the rest of the run
    // down with it. The cost of 'ignore' is that workflow.success stays
    // true no matter how many predictions were dropped.
    //
    // So count rows rather than trust the exit status. Both files carry a
    // header row, hence the -1. A null return means the file is absent,
    // which for the rankings file IS the diagnosis.
    // -------------------------------------------------------------------
    def dataRows = { pathStr ->
        def f = file(pathStr)
        f.exists() ? Math.max(0, f.readLines().count { it.trim() } - 1) : null
    }

    def n_pairs  = dataRows("${params.outdir}/pairs/pairs.tsv")
    def n_ranked = dataRows("${params.outdir}/rankings_${params.predictor}.tsv")

    // Under --msa_only there is no rankings file by design, so reconcile
    // against the thing that run was actually for: one cached MSA output
    // per protein the pair list names.
    def msaReport = {
        def pf = file("${params.outdir}/pairs/pairs.tsv")
        if (!pf.exists()) return null
        def lines = pf.readLines().findAll { it.trim() }
        if (lines.size() < 2) return [0, 0]
        def hdr = lines[0].split('\t') as List
        def bi = hdr.indexOf('bait_id'), ci = hdr.indexOf('candidate_id')
        if (bi < 0 || ci < 0) return null
        def ids = [] as Set
        lines.drop(1).each { def f = it.split('\t'); ids << f[bi]; ids << f[ci] }
        def af3 = params.predictor.toLowerCase() == 'alphafold3'
        def dir = af3 ? params.alignments_dir_af3 : params.alignments_dir_boltz2
        def present = ids.count { id ->
            file(af3 ? "${dir}/a3m/${id}.a3m" : "${dir}/csv/${id}.csv").exists()
        }
        return [present, ids.size()]
    }

    def n_ignored = workflow.stats.ignoredCount

    def status
    def msa = params.msa_only ? msaReport() : null
    if (!workflow.success)
        // A real pipeline failure (not an ignored task): Nextflow has
        // already printed the cause above.
        status = 'FAILED'
    else if (params.msa_only && msa == null)
        status = 'SUCCESS (--msa_only; no pair list found to reconcile against)'
    else if (params.msa_only && msa[0] < msa[1])
        status = "INCOMPLETE - ${msa[0]} of ${msa[1]} proteins have MSAs (--msa_only; no predictions run)"
    else if (params.msa_only)
        status = "SUCCESS - ${msa[0]} of ${msa[1]} proteins have MSAs (--msa_only; no predictions run)"
    else if (n_pairs == null)
        // No pairs.tsv published, so there is nothing to reconcile against
        // (e.g. GENERATE_PAIRS itself was the thing that died).
        status = 'SUCCESS (no pair list found to reconcile against)'
    else if (n_ranked == null)
        status = "FAILED - no rankings written (0 of ${n_pairs} pairs ranked)"
    else if (n_ranked < n_pairs)
        status = "INCOMPLETE - ${n_ranked} of ${n_pairs} pairs ranked, ${n_pairs - n_ranked} dropped"
    else
        status = "SUCCESS - ${n_ranked} of ${n_pairs} pairs ranked"

    // Assembled line by line rather than as one GString so the optional
    // "Dropped" line can be omitted without leaving a blank row or
    // disturbing stripIndent()'s common-prefix calculation.
    def lines = []
    lines << '========================================================================='
    lines << 'OmegaPulldown complete'
    lines << '========================================================================='
    lines << "Predictor  : ${params.predictor}"
    lines << "Status     : ${status}"
    lines << "Duration   : ${workflow.duration}"
    if (n_ignored > 0)
        lines << "Dropped    : ${n_ignored} task(s) failed after retries and were ignored"
    lines << "Output     : ${params.outdir}"
    if (params.msa_only) {
        // Name the cache, not the structures/rankings this run never wrote.
        lines << "Alignments : " + (params.predictor.toLowerCase() == 'alphafold3'
                                    ? params.alignments_dir_af3
                                    : params.alignments_dir_boltz2)
    }
    else {
        lines << "Structures : ${params.outdir}/structures_${params.predictor}/"
        lines << "Rankings   : ${params.outdir}/rankings_${params.predictor}.tsv"
    }
    lines << '========================================================================='

    log.info '\n' + lines.collect { '    ' + it }.join('\n') + '\n'
}
