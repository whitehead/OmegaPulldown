/*
 * AF3-native MSA search (predictor=alphafold3 only).
 *
 * Replicates AlphaFold 3's own genetic search using mmseqs2 as a validated
 * drop-in for jackhmmer/HHblits (see the AlphaFast paper, Perry,Kim and Romero
 * 2026), against the exact database versions AF3 was trained/built around.
 *
 * Two independent outputs per unique protein, matching AF3-native's own
 * architecture (confirmed from RomeroLab/alphafast's fork of
 * alphafold3/data/pipeline.py):
 *   - Unpaired MSA: UniRef90 + MGnify + Small BFD, merged with
 *     exact-duplicate-sequence removal (MSA_SEARCH_DB + MERGE_AF3_NATIVE_MSA).
 *   - Paired MSA: UniProt only, kept completely separate, full original
 *     headers preserved so AlphaFold 3's own species-pairing code
 *     (alphafold3/data/msa_identifiers.py) can extract per-hit species IDs
 *     and do bait/candidate cross-chain matching itself at inference time
 *     (MSA_SEARCH_UNIPROT_PAIRED). No pairing logic is implemented here:
 *     AF3 does it internally given a correctly-headered pairedMsaPath.
 *
 * db is a map: [name: <string>, path: <mmseqs2 db path>, max_seqs: <int>]
 * (see params.af3_native_dbs in nextflow.config). One invocation per
 * (batch-of-proteins, database): each database searched independently
 * against the raw queries, unlike ColabFold's cascaded profile search; no
 * profile is carried between databases, matching AF3-native's own
 * methodology.
 *
 * Batched across proteins: the expensive part of this process is loading
 * the large mmseqs2 index (hundreds of GB for UniRef90/MGnify), not the
 * search itself for one short query. The multi-query result is split back
 * into per-protein files with `mmseqs unpackdb`.
 */
process MSA_SEARCH_DB {
    tag "batch_of_${batch.size()}:${db.name}"
    label 'process_mmseqs_large'

    input:
    tuple val(batch), val(db)

    output:
    path("*.${db.name}.a3m"), emit: a3m_files

    script:
    // batch is a list of [prot_id, sequence] pairs (see main.nf's
    // split_ch.compute.collate(batch_size) -> combine(db_ch)).
    def fasta_lines = batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')
    // --split-memory-limit derived from the task's actual resolved memory,
    // not hardcoded, so it can't drift out of sync with it. 10GB headroom
    // for non-index overhead (per-thread buffers, I/O, result buffers) that
    // --split alone doesn't account for.
    def split_mem = "${task.memory.toGiga() - 10}G"
    """
    # Write query FASTA: one record per protein in this batch
    cat << 'EOF' > query.fasta
${fasta_lines}
EOF

    # Create query database (one multi-record DB for the whole batch)
    mmseqs createdb query.fasta queryDB

    # Search against this single AF3-native database.
    mmseqs search queryDB ${db.path} resultDB tmp \
        --split-memory-limit ${split_mem} \
        --split ${db.split} \
        --split-mode 0 \
        -a \
        -s ${params.mmseqs_sensitivity} \
        -e ${params.af3_native_evalue} \
        --max-seqs ${db.max_seqs} \
        --threads ${task.cpus}

    # Convert to a3m MSA format (one packed multi-query result)
    mmseqs result2msa queryDB ${db.path} resultDB packed.a3m \
        --msa-format-mode 6

    # Split the packed result back into one file per protein.
    # --unpack-name-mode 1 names each file by "accession, through .lookup"
    cp queryDB.lookup packed.a3m.lookup
    mmseqs unpackdb packed.a3m . --unpack-name-mode 1 --unpack-suffix ".${db.name}.a3m"

    # Clean up tmp
    rm -rf tmp queryDB* resultDB* packed.a3m packed.a3m.lookup
    """

    stub:
    // No mmseqs invocation: avoids real compute against production DBs
    // under -stub-run.
    def stub_lines = batch.collect { prot_id, sequence ->
        "echo '>${prot_id}' > ${prot_id}.${db.name}.a3m\necho '${sequence}' >> ${prot_id}.${db.name}.a3m"
    }.join('\n')
    """
    ${stub_lines}
    """
}

/*
 * Merge the three per-database a3m files (UniRef90, MGnify, Small BFD) for
 * one protein into a single unpaired MSA, with exact-duplicate-sequence
 * removal (see merge_af3_native_msa.py for the AF3-native dedup semantics).
 * UniProt is NOT one of the inputs, see MSA_SEARCH_UNIPROT_PAIRED.
 */
process MERGE_AF3_NATIVE_MSA {
    tag "${prot_id}"
    label 'process_low'

    // params.alignments_dir_af3, NOT params.outdir. This directory IS the
    // per-protein MSA cache: main.nf's .exists() branch looks each output
    // up under params.alignments_dir_af3, so publishing anywhere else means
    // writing where nothing reads.
    //
    publishDir "${params.alignments_dir_af3}/a3m", mode: 'copy'

    input:
    tuple val(prot_id), path(a3m_files)

    output:
    tuple val(prot_id), path("${prot_id}.a3m"), emit: a3m

    script:
    """
    merge_af3_native_msa.py ${a3m_files} --output ${prot_id}.a3m
    """

    stub:
    // Pure-Python, local files only, no heavy/external dependency: safe
    // and fast to run for real even under -stub-run.
    """
    merge_af3_native_msa.py ${a3m_files} --output ${prot_id}.a3m
    """
}

/*
 * UniProt search for the paired MSA, batched across proteins (see
 * MSA_SEARCH_DB above for the batching rationale, identical mechanism
 * here). Uses --msa-format-mode 2 (aligned FASTA, full original headers
 * preserved) instead of mode 6 (used by MSA_SEARCH_DB above): mode 6
 * discards the original header text in favor of a synthetic ID+stats
 * header, but AF3's own species-pairing code needs the raw UniProt header
 * (sp|ACCESSION|GENE_SPECIES / tr|ACCESSION|GENE_SPECIES) to extract a
 * per-hit species identifier. e-value/max-seqs match AF3-native's own
 * UniProt defaults (max_sequences=50,000, the deepest of the four DBs).
 */
process MSA_SEARCH_UNIPROT_PAIRED {
    tag "batch_of_${batch.size()}:uniprot_paired"
    label 'process_medium'

    // params.alignments_dir_af3, not params.outdir: see the identical
    // note on MERGE_AF3_NATIVE_MSA above. The paired a3m is cached and
    // looked up independently of the unpaired one (main.nf:505), so it has
    // to track the same parameter or half the cache goes missing.
    publishDir "${params.alignments_dir_af3}/paired_a3m", mode: 'copy'

    input:
    val(batch)

    output:
    // Named "<prot_id>.paired.a3m", not "<prot_id>.a3m": these files get
    // staged into the same MAKE_AF3_JSON task directory as the unpaired
    // "<prot_id>.a3m" for the same protein.
    path("*.paired.a3m"), emit: a3m_files

    script:
    def fasta_lines = batch.collect { prot_id, sequence -> ">${prot_id}\n${sequence}" }.join('\n')
    def split_mem = "${task.memory.toGiga() - 10}G"
    """
    # Write query FASTA, one record per protein in this batch
    cat << 'EOF' > query.fasta
${fasta_lines}
EOF

    # Create query database (one multi-record DB for the whole batch)
    mmseqs createdb query.fasta queryDB

    # --split + --split-mode 0 explicit for the same reason as MSA_SEARCH_DB
    # above: --split-mode defaults to 2 (auto), which ignores an explicit
    # --split entirely.
    mmseqs search queryDB ${params.af3_native_uniprot_db.path} resultDB tmp \
        --split-memory-limit ${split_mem} \
        --split ${params.af3_native_uniprot_db.split} \
        --split-mode 0 \
        -a \
        -s ${params.mmseqs_sensitivity} \
        -e ${params.af3_native_evalue} \
        --max-seqs ${params.af3_native_uniprot_db.max_seqs} \
        --threads ${task.cpus}

    # Aligned FASTA (mode 2), not A3M (mode 6): preserves full headers.
    # One packed multi-query result for the whole batch.
    mmseqs result2msa queryDB ${params.af3_native_uniprot_db.path} resultDB packed_aligned.fasta \
        --msa-format-mode 2

    # Split the packed result back into one aligned-FASTA file per protein,
    # named by prot_id (--unpack-name-mode 1, see MSA_SEARCH_DB above for why
    # mode 0 is wrong here: it names by internal DB key, not prot_id; and
    # why the lookup file must be copied under packed_aligned.fasta's own
    # name first: result2msa's output doesn't inherit queryDB.lookup).
    cp queryDB.lookup packed_aligned.fasta.lookup
    mmseqs unpackdb packed_aligned.fasta . --unpack-name-mode 1 --unpack-suffix ".aligned.fasta"

    # Per-protein a3m conversion -- cheap pure-Python step, done in a loop
    # here rather than fanning back out to one task per protein.
    for f in *.aligned.fasta; do
        prot_id="\${f%.aligned.fasta}"
        aligned_fasta_to_a3m.py "\$f" --output "\${prot_id}.paired.a3m"
    done

    # Clean up tmp
    rm -rf tmp queryDB* resultDB* packed_aligned.fasta packed_aligned.fasta.lookup *.aligned.fasta
    """

    stub:
    // No mmseqs invocation: avoids real compute against production DBs
    // under -stub-run. One query-only a3m (no hits) per protein in the
    // batch, a valid minimal MSA.
    def stub_lines = batch.collect { prot_id, sequence ->
        "echo '>${prot_id}' > ${prot_id}.paired.a3m\necho '${sequence}' >> ${prot_id}.paired.a3m"
    }.join('\n')
    """
    ${stub_lines}
    """
}
