#!/usr/bin/env python3
"""Retained spatial effects and true-donor LOO; no new p values or FDR.

Same-bin: log2 of the unweighted median section g, within donor or after
omitting all sections from one donor. Original family/headline labels are
historical labels only. Five LOO folds summarize influence, not confidence.
P16 distance summary: symmetric directional log2(g) per ring, ring-width
weighted within section, then an unweighted median across sections. The
same-bin thresholds are not applied to the three distance summaries.
The original P1-P36/P13-P54 regional models and intervals are untouched.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__import__("os").environ["NMF_WORK_ROOT"])
ARCHIVE = ROOT / 'inputs/cortex_nmf_program/archived'
CROSS = ARCHIVE / 'results/crossregion_v1'
PCDIR = CROSS / 'markcorr_v2/per_chip'
DONOR_TABLE = ARCHIVE / 'revisions/2026-08-02_v49_professor_review/human_validation/HUMAN_VALIDATION_spatial_section_by_layer_aggregates_all54.tsv'
OUTPUT = ROOT / 'analysis/02_true_donor_spatial_support'
THRESHOLD = 0.32
STABLE_FRACTION = 0.90
P16_PAIRS = [('P16', 'P36'), ('P16', 'P45'), ('P16', 'P54')]


def raw_id(label):
    return int(str(label).removeprefix('program_'))


def log_median(g):
    median = float(np.median(g))
    # Same-bin source definition: nonpositive median has no finite log2 effect.
    return float(np.log2(median)) if median > 0 else float('nan')


def sign_fraction(values, reference):
    return float(np.mean(np.sign(values) == np.sign(reference))) if reference != 0 and np.isfinite(reference) else float('nan')


def main():
    mapping = pd.read_csv(CROSS / 'program_renumber_map.tsv', sep='\t')
    mapping = mapping.loc[mapping.status.eq('kept')].copy()
    mapping['new_int'] = mapping.new_P.astype(str).str.removeprefix('P').astype(int)
    mapping = mapping.sort_values('new_int')
    old_to_new = dict(zip(mapping.old_P.astype(int), 'P' + mapping.new_int.astype(str)))
    new_to_old = {new: old for old, new in old_to_new.items()}
    metadata = pd.read_csv(DONOR_TABLE, sep='\t', usecols=['chip', 'donor', 'region']).drop_duplicates()
    if metadata.chip.duplicated().any() or metadata[['chip', 'donor', 'region']].isna().any().any():
        raise ValueError('Ambiguous or missing chip/donor/region key in existing source')
    metadata = metadata.set_index('chip')
    donors = sorted(metadata.donor.unique())
    history = pd.read_csv(ROOT / 'tables/TableS6_between_chip_colocalization.tsv', sep='\t')
    historical_flags = {}
    for row in history.itertuples():
        mode = str(row.mode)
        if mode == 'progprog':
            a, b = old_to_new[int(row.A_cnmf_id)], old_to_new[int(row.B_cnmf_id)]
            a, b = sorted((a, b), key=lambda x: int(x[1:]))
        elif mode == 'cellprog':
            a, b = str(row.A_label), old_to_new[int(row.B_cnmf_id)]
        else:
            continue
        historical_flags[(mode, a, b)] = str(row.is_headline).lower() == 'true'

    donor_rows, loo_rows, ring_rows = [], [], []
    p16_sections = {pair: [] for pair in P16_PAIRS}
    for mode in ['cellprog', 'progprog']:
        chips, same_bin = [], []
        pairs = None
        reference_axes = None
        for path in sorted(PCDIR.glob(mode + '_*.npz')):
            with np.load(path, allow_pickle=True) as data:
                chip = str(data['chip'])
                meta = metadata.loc[chip]
                a_names = [str(x) for x in data['A_names']]
                b_names = [str(x) for x in data['B_names']]
                edges = np.asarray(data['ring_edges_um'], dtype=float)
                if reference_axes is None:
                    reference_axes = (a_names, b_names, edges.tolist())
                    b_index = {raw_id(name): i for i, name in enumerate(b_names)}
                    if mode == 'cellprog':
                        pairs = [(ai, b_index[old], a, new, '', old)
                                 for ai, a in enumerate(a_names) for old, new in old_to_new.items()]
                    else:
                        a_index = {raw_id(name): i for i, name in enumerate(a_names)}
                        retained = list(old_to_new)
                        pairs = [(a_index[oa], b_index[ob], old_to_new[oa], old_to_new[ob], oa, ob)
                                 for i, oa in enumerate(retained) for ob in retained[i + 1:]]
                elif reference_axes != (a_names, b_names, edges.tolist()):
                    raise ValueError('Source tensor axes differ between sections')
                if edges.tolist() != [0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500]:
                    raise ValueError('Unexpected source ring boundaries')
                g = np.nan_to_num(data['g'].astype(float), nan=1.0, posinf=1.0, neginf=1.0)
                chips.append(chip)
                same_bin.append([g[ai, bi, 0] for ai, bi, *_ in pairs])
                if mode == 'progprog':
                    logg = np.log2(np.clip(g, 1e-9, None))
                    a_index = {raw_id(name): i for i, name in enumerate(a_names)}
                    b_index = {raw_id(name): i for i, name in enumerate(b_names)}
                    widths = np.diff(edges)
                    for pair in P16_PAIRS:
                        left, right = pair
                        oa, ob = new_to_old[left], new_to_old[right]
                        lr = logg[a_index[oa], b_index[ob], :]
                        rl = logg[a_index[ob], b_index[oa], :]
                        symmetric = (lr + rl) / 2.0
                        # Same left-to-right summation order as distance_weighted_mean.
                        integrated = sum(float(v) * float(w) for v, w in zip(symmetric, widths)) / sum(float(w) for w in widths)
                        p16_sections[pair].append((chip, str(meta.donor), str(meta.region), integrated))
                        for r in range(len(widths)):
                            ring_rows.append(dict(pair=f'{left}-{right}', left_program=left, right_program=right,
                                left_component=oa, right_component=ob, chip=chip, donor=str(meta.donor), region=str(meta.region),
                                ring_start_um=edges[r], ring_end_um=edges[r + 1], ring_width_um=widths[r],
                                log2g_left_to_right=lr[r], log2g_right_to_left=rl[r], symmetric_log2g=symmetric[r],
                                section_relation_0_500=integrated))
        chips = np.asarray(chips)
        g_matrix = np.asarray(same_bin, dtype=float)
        chip_donors = metadata.loc[chips, 'donor'].to_numpy()
        for j, (_, _, a, b, oa, ob) in enumerate(pairs):
            values = g_matrix[:, j]
            full = log_median(values)
            folds = np.asarray([log_median(values[chip_donors != donor]) for donor in donors])
            retained = np.isfinite(folds) & (np.abs(folds) > THRESHOLD) & (np.sign(folds) == np.sign(full))
            fraction = float(np.mean(retained))
            common = dict(mode=mode, A_label=a, B_label=b, A_component=oa, B_component=ob,
                historical_is_headline=historical_flags[(mode, a, b)],
                n_sections_all=len(chips), n_donors_all=len(donors), all_sections_log2_median_g=full)
            for k, donor in enumerate(donors):
                inside = chip_donors == donor
                donor_rows.append(dict(**common, donor=donor, n_sections_donor=int(inside.sum()),
                    donor_log2_median_g=log_median(values[inside])))
                loo_rows.append(dict(**common, omitted_donor=donor,
                    n_sections_omitted=int(inside.sum()), n_sections_remaining=int((~inside).sum()),
                    n_donors_remaining=len(donors) - 1, loo_log2_median_g=folds[k],
                    loo_same_direction=bool(np.isfinite(folds[k]) and np.sign(folds[k]) == np.sign(full)),
                    loo_retains_original_effect_rule=bool(retained[k]),
                    loo_min_log2_median_g=float(np.min(folds)), loo_max_log2_median_g=float(np.max(folds)),
                    frac_folds_retaining_original_effect_rule=fraction,
                    stable_90=bool(fraction >= STABLE_FRACTION)))

    p16_donor_rows, p16_loo_rows = [], []
    for pair, records in p16_sections.items():
        frame = pd.DataFrame(records, columns=['chip', 'donor', 'region', 'section_relation_0_500'])
        values = frame.section_relation_0_500.to_numpy()
        full = float(np.median(values))
        folds = np.asarray([np.median(frame.loc[frame.donor.ne(d), 'section_relation_0_500']) for d in donors])
        common = dict(pair='-'.join(pair), left_program=pair[0], right_program=pair[1],
            n_sections_all=len(frame), n_donors_all=len(donors),
            all_sections_median_relation_0_500=full,
            all_sections_same_sign_fraction=sign_fraction(values, full))
        for k, donor in enumerate(donors):
            inside = frame.donor.eq(donor)
            donor_values = frame.loc[inside, 'section_relation_0_500'].to_numpy()
            donor_effect = float(np.median(donor_values))
            p16_donor_rows.append(dict(**common, donor=donor, n_sections_donor=int(inside.sum()),
                donor_median_relation_0_500=donor_effect,
                donor_sections_same_sign_as_all_fraction=sign_fraction(donor_values, full)))
            p16_loo_rows.append(dict(**common, omitted_donor=donor,
                n_sections_omitted=int(inside.sum()), n_sections_remaining=int((~inside).sum()),
                n_donors_remaining=len(donors) - 1, loo_median_relation_0_500=float(folds[k]),
                loo_same_direction=bool(np.sign(folds[k]) == np.sign(full)),
                loo_min_median_relation_0_500=float(np.min(folds)),
                loo_max_median_relation_0_500=float(np.max(folds)),
                frac_folds_same_direction=sign_fraction(folds, full)))

    OUTPUT.mkdir(parents=True, exist_ok=True)
    outputs = {
        'same_bin_per_donor_effects.tsv': pd.DataFrame(donor_rows),
        'same_bin_true_donor_loo.tsv': pd.DataFrame(loo_rows),
        'P16_distance_0_500_per_donor_effects.tsv': pd.DataFrame(p16_donor_rows),
        'P16_distance_0_500_true_donor_loo.tsv': pd.DataFrame(p16_loo_rows),
        'P16_section_ring_effects.tsv': pd.DataFrame(ring_rows),
    }
    for name, table in outputs.items():
        table.to_csv(OUTPUT / name, sep='\t', index=False)
        print(f'RETAINED {OUTPUT / name}: {len(table)} rows', flush=True)
    print('True donor section counts:', metadata.groupby('donor').size().to_dict(), flush=True)
    same = outputs['same_bin_true_donor_loo.tsv'].drop_duplicates(['mode', 'A_label', 'B_label'])
    for mode, frame in same.groupby('mode'):
        old = frame.loc[frame.historical_is_headline]
        print(f'{mode}: {len(frame)} pairs; historical headline {len(old)}; true-donor stable_90 within historical headline {int(old.stable_90.sum())}', flush=True)
    for pair in P16_PAIRS:
        name = '-'.join(pair)
        dt = outputs['P16_distance_0_500_per_donor_effects.tsv']
        dt = dt.loc[dt.pair.eq(name)]
        lt = outputs['P16_distance_0_500_true_donor_loo.tsv']
        lt = lt.loc[lt.pair.eq(name)].iloc[0]
        print(name, 'all-section effect', lt.all_sections_median_relation_0_500,
              'donor effects', dict(zip(dt.donor, dt.donor_median_relation_0_500)),
              'LOO range', lt.loo_min_median_relation_0_500, lt.loo_max_median_relation_0_500,
              'LOO direction fraction', lt.frac_folds_same_direction, flush=True)


if __name__ == '__main__':
    main()
