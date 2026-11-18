/*
 * Parse AlphaFold 3 outputs and rank predictions by interface quality.
 *
 * AF3 output structure per pair:
 *   <pair_id>/
 *     <pair_id>_summary_confidences.json   (top-ranked model summary)
 *     <pair_id>_ranking_scores.csv         (per-sample ranking scores)
 *     <pair_id>_model.cif                  (top-ranked structure)
 *     seed-<N>_sample-<M>/
 *       <pair_id>_seed-<N>_sample-<M>_model.cif
 *       <pair_id>_seed-<N>_sample-<M>_summary_confidences.json
 *
 * Reuses RANK_PREDICTIONS from rank.nf — the output TSV format is identical.
 */

process PARSE_METRICS_AF3 {
    tag "${pair_id}"
    label 'process_low'

    input:
    tuple val(pair_id), path(pred_files)

    output:
    path "${pair_id}.metrics.tsv", emit: metrics

    script:
    """
    #!/usr/bin/env python3
    import json, glob, math, re, sys, os, shutil
    import numpy as np

    pair_id = "${pair_id}"

    # pDockQ comes from bin/compute_pdockq.py.
    _p = shutil.which("compute_pdockq.py")
    if _p is None:
        sys.exit("compute_pdockq.py not found on PATH -- is bin/ present?")
    sys.path.insert(0, os.path.dirname(_p))
    from compute_pdockq import (compute_pdockq, parse_cif_atoms,
                                get_cb_atoms, INTERFACE_CUTOFF)

    def compute_plddts(cif_path):
        \"\"\"
        True pLDDT metrics for an AF3 model.

        AF3's summary_confidences.json carries no pLDDT at all, but AF3
        writes per-atom pLDDT (0-100) into the CIF B-factor column, so both
        metrics are derived from the structure:

          complex_plddt  : mean pLDDT over every atom in the model.
          complex_iplddt : mean pLDDT over every atom of the interface
                           residues. A residue is at the interface when its
                           CB (CA for glycine) lies within INTERFACE_CUTOFF
                           of a CB on the other chain: the same definition that
                           compute_pdockq uses, so the two metrics agree on
                           what "interface" means.

        Do not compare either number against Boltz-2's complex_plddt,
        which is a TOKEN (per-residue) mean.

        complex_iplddt is deliberately not the same number as
        avg_interface_plddt: that averages the single representative CB
        per interface residue (residue-weighted), this averages all atoms of
        those residues (atom-weighted), so side-chain confidence is
        reflected. Same interface, two resolutions, two weightings.

        Returns (complex_plddt, complex_iplddt), using None for a metric that
        is undefined: a monomer, or a model with no cross-chain contact.
        None is not 0.0: "no interface" must not read as "an interface that
        scores zero", which is the distinction compute_pdockq already keeps
        for pdockq.
        \"\"\"
        atoms = parse_cif_atoms(cif_path)
        if not atoms:
            return None, None

        complex_plddt = float(np.mean([a["bfactor"] for a in atoms]))

        cb = get_cb_atoms(atoms)
        chains = sorted(set(k[0] for k in cb.keys()))
        if len(chains) < 2:
            return complex_plddt, None

        res_a = [k for k in cb if k[0] == chains[0]]
        res_b = [k for k in cb if k[0] == chains[1]]
        if not res_a or not res_b:
            return complex_plddt, None

        xyz_a = np.array([cb[k][:3] for k in res_a])
        xyz_b = np.array([cb[k][:3] for k in res_b])
        dist = np.sqrt(((xyz_a[:, None, :] - xyz_b[None, :, :]) ** 2).sum(axis=2))
        mask = dist < INTERFACE_CUTOFF
        if not mask.any():
            return complex_plddt, None

        iface = set()
        for idx in np.where(mask.any(axis=1))[0]:
            iface.add(res_a[idx])
        for idx in np.where(mask.any(axis=0))[0]:
            iface.add(res_b[idx])

        vals = [a["bfactor"] for a in atoms
                if (a["chain_id"], a["res_seq"]) in iface]
        if not vals:
            return complex_plddt, None
        return complex_plddt, float(np.mean(vals))

    # =====================================================================
    # Collect per-sample metrics from AF3 output
    # =====================================================================
    # AF3 metrics mapping to Boltz-2 column names for unified ranking:
    #   AF3 ranking_score -> af3_ranking_score  (NOT Boltz-2 formula!)
    #   AF3 ptm           -> ptm
    #   AF3 iptm          -> iptm  (also used as protein_iptm)
    #   AF3 CIF B-factors -> af3_complex_plddt, af3_complex_iplddt
    #   AF3 has_clash     -> not read
    #
    # The two rankings tables share the COMPARABLE part of their schema and
    # nothing else. The column names say which is which: a bare name
    # is shared and meant to be read across both tables, while an
    # af3_ / boltz2_ prefix marks a number that is only meaningful within its
    # own table.
    #
    #       af3_ranking_score       = 0.8*iptm + 0.2*ptm
    #                                 + 0.5*fraction_disordered
    #                                 - 100*has_clash
    #       boltz2_confidence_score = 0.8*complex_plddt + 0.2*iptm   (to 1e-8)
    #
    # Weighting conventions inside the AF3 row, since they are not uniform:
    # af3_complex_plddt and af3_complex_iplddt are atom-weighted,
    # avg_interface_plddt is residue-weighted (one representative CB per
    # interface residue).
    conf_keys = ['af3_ranking_score','ptm','iptm','protein_iptm',
                 'af3_complex_plddt','af3_complex_iplddt']
    struct_keys = ['pdockq','n_interface_contacts','avg_interface_plddt']
    all_keys = conf_keys + struct_keys

    sample_vals = {k: [] for k in all_keys}

    # Find per-sample confidence JSONs
    sample_conf_files = sorted(glob.glob("*_sample-*_summary_confidences.json"))
    sample_cif_files  = sorted(glob.glob("*_sample-*_model.cif"))

    # If no per-sample files, fall back to the top-level summary
    if not sample_conf_files:
        sample_conf_files = sorted(glob.glob("*_summary_confidences.json"))
        # Exclude any that contain "sample-" (shouldn't happen, but be safe)
        sample_conf_files = [f for f in sample_conf_files if "sample-" not in f]
    if not sample_cif_files:
        sample_cif_files = sorted(glob.glob("*_model.cif"))
        sample_cif_files = [f for f in sample_cif_files if "sample-" not in f]

    n_samples = max(len(sample_conf_files), len(sample_cif_files), 1)

    for i in range(n_samples):
        # --- Confidence metrics ---
        if i < len(sample_conf_files):
            with open(sample_conf_files[i]) as f:
                data = json.load(f)

            # AF3 -> unified column mapping
            if 'ranking_score' in data:
                sample_vals['af3_ranking_score'].append(float(data['ranking_score']))
            if 'ptm' in data:
                sample_vals['ptm'].append(float(data['ptm']))
            if 'iptm' in data:
                sample_vals['iptm'].append(float(data['iptm']))
                # AF3 iptm is the inter-chain metric; map to protein_iptm for ranking
                sample_vals['protein_iptm'].append(float(data['iptm']))

        # --- Structural metrics (pDockQ, pLDDT) ---
        if i < len(sample_cif_files):
            pdockq, nc, aip = compute_pdockq(sample_cif_files[i])
            if pdockq is not None:
                sample_vals['pdockq'].append(pdockq)
            sample_vals['n_interface_contacts'].append(float(nc))
            sample_vals['avg_interface_plddt'].append(aip)

            cplddt, iplddt = compute_plddts(sample_cif_files[i])
            if cplddt is not None:
                sample_vals['af3_complex_plddt'].append(cplddt)
            if iplddt is not None:
                sample_vals['af3_complex_iplddt'].append(iplddt)

    if not sample_conf_files and not sample_cif_files:
        print(f"WARNING: no confidence JSON or CIF for {pair_id}", file=sys.stderr)

    # =====================================================================
    # Compute mean and std, write TSV (same format as Boltz-2 metrics)
    # =====================================================================
    out_cols = ['pair_id', 'n_samples']
    for k in all_keys:
        out_cols.append(f"{k}_mean")
        out_cols.append(f"{k}_std")

    out_vals = [pair_id, str(n_samples)]
    for k in all_keys:
        arr = sample_vals[k]
        if len(arr) == 0:
            out_vals.extend(["NA", "NA"])
        else:
            a = np.array(arr)
            out_vals.append(f"{np.mean(a):.4f}")
            out_vals.append(f"{np.std(a, ddof=0):.4f}" if len(a) > 1 else "NA")

    with open("${pair_id}.metrics.tsv", "w") as out:
        out.write("\\t".join(out_cols) + "\\n")
        out.write("\\t".join(out_vals) + "\\n")

    print(f"{pair_id}: {n_samples} samples, "
          f"iptm = {np.mean(sample_vals['iptm']):.4f} "
          f"+/- {np.std(sample_vals['iptm'], ddof=0):.4f}"
          if sample_vals['iptm'] else f"{pair_id}: no confidence data",
          file=sys.stderr)
    """

    stub:
    // Minimal placeholder, not a duplicate of the embedded Python parser
    // above -- just satisfies the declared output for wiring checks.
    """
    printf "pair_id\\tn_samples\\taf3_ranking_score_mean\\taf3_ranking_score_std\\n" > ${pair_id}.metrics.tsv
    printf "${pair_id}\\t1\\t0.5\\tNA\\n" >> ${pair_id}.metrics.tsv
    """
}
