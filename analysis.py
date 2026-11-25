"""Evaluate three precomputed linear B-cell scores against ovalbumin assays."""
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

from common import DATA, RESULTS, plot_style, save_summary

METHODS = ("BepiPred-2.0", "Parker", "Emini")
CONTEXT = ("pubmed_id", "assay_names", "host_organism_name", "antibody_isotype")


def score_tracks():
    reference = "".join((DATA / "ova.fasta").read_text().splitlines()[1:])
    frame = pd.read_csv(DATA / "external_scores.csv")
    if "".join(frame.aa) != reference or list(frame.position) != list(range(1, len(reference) + 1)):
        raise ValueError("Predictor coordinates differ from UniProt")
    for method in METHODS:
        frame[method + "_percentile"] = frame[method].rank(method="average", pct=True)
    return reference, frame


def audit(assays, tracks, reference):
    assays = assays[assays.assay_names.isin([
        "qualitative binding|microarray", "qualitative binding|western blot"
    ])].copy()
    assays["start"] = assays.linear_sequence.map(reference.find)
    assays["matches"] = assays.linear_sequence.map(reference.count)
    mapped = assays[assays.matches.eq(1)]
    mapped = mapped.assign(outcome=mapped.qualitative_measure.str.startswith("Positive").astype(int))
    rows, conflicts, missing_scores = [], 0, 0
    for key, group in mapped.groupby(list(CONTEXT), dropna=False):
        tests = group.groupby("linear_sequence").outcome.agg(["min", "max", "size"])
        valid = tests[tests["min"].eq(tests["max"])]
        conflicts += len(tests) - len(valid)
        for sequence, row in valid.iterrows():
            first = reference.index(sequence)
            block = tracks.iloc[first:first+len(sequence)]
            if block[list(METHODS)].isna().any().any():
                missing_scores += 1
                continue
            rows.append({
                **dict(zip(CONTEXT, key)), "sequence": sequence, "start": first + 1,
                "end": first + len(sequence), "positive": int(row["min"]), "replicate_rows": int(row["size"]),
                **{name: float(block[name].mean()) for name in METHODS},
            })
    return pd.DataFrame(rows), len(assays) - len(mapped), conflicts, missing_scores


def evaluate(table):
    rows = []
    for identifiers, group in table.groupby(list(CONTEXT), dropna=False):
        if group.positive.nunique() != 2 or min(group.positive.value_counts()) < 3:
            continue
        for method in METHODS:
            rows.append({
                **dict(zip(CONTEXT, identifiers)), "method": method, "peptides": len(group),
                "positive": int(group.positive.sum()), "prevalence_ap_baseline": group.positive.mean(),
                "roc_auc": roc_auc_score(group.positive, group[method]),
                "average_precision": average_precision_score(group.positive, group[method]),
            })
    return pd.DataFrame(rows)


def main():
    reference, tracks = score_tracks()
    assays = pd.read_csv(DATA / "assays.csv").fillna("")
    if set(assays.parent_source_antigen_iri) != {"UNIPROT:P01012"}:
        raise ValueError("Unexpected antigen in assay table")
    table, missing, conflicts, missing_scores = audit(assays, tracks, reference)
    metrics = evaluate(table)
    if metrics.empty:
        raise ValueError("No comparable assay contexts")
    correlations = pd.DataFrame([
        {"a": a, "b": b, "spearman_rho": spearmanr(tracks[a], tracks[b], nan_policy="omit").statistic}
        for a, b in combinations(METHODS, 2)
    ])
    RESULTS.mkdir(exist_ok=True)
    tracks.to_csv(RESULTS / "region_scores.csv", index=False, float_format="%.6f")
    table.to_csv(RESULTS / "assayed_regions.csv", index=False, float_format="%.6f")
    metrics.to_csv(RESULTS / "metrics.csv", index=False, float_format="%.6f")
    correlations.to_csv(RESULTS / "model_agreement.csv", index=False, float_format="%.6f")
    summary = metrics.groupby("method")[["roc_auc", "average_precision"]].mean()
    summary.to_csv(RESULTS / "macro_average.csv", float_format="%.6f")

    plt = plot_style()
    fig, (ax, bx) = plt.subplots(
        2, 1, figsize=(10, 6.5), gridspec_kw={"height_ratios": [2.0, 1.2]}, constrained_layout=True
    )
    for name, color in zip(METHODS, ["#1f766e", "#c17049", "#4e69a2"]):
        ax.plot(tracks.position, tracks[name + "_percentile"], label=name, color=color, alpha=0.82)
    ax.set(xlabel="UniProt P01012 position", ylabel="Within-protein percentile",
           title="Predictor profiles across chicken ovalbumin", ylim=(0, 1.02))
    ax.legend(frameon=False, ncol=3)
    summary.plot.barh(ax=bx, color=["#1f766e", "#4e69a2"])
    bx.axvline(0.5, linestyle=":", color="#7c8992")
    bx.set(xlabel="Macro average across assay contexts", ylabel="", xlim=(0, 1.02))
    bx.legend(frameon=False, fontsize=8)
    fig.savefig(RESULTS / "predictor_comparison.png")
    plt.close(fig)
    save_summary({
        "closed_reference": "UniProt P01012", "source_assay_records": len(assays),
        "retained_peptide_context_pairs": len(table), "unmapped_excluded": missing,
        "conflicting_context_pairs_excluded": conflicts,
        "peptide_contexts_missing_terminal_scores": missing_scores,
        "scored_contexts": len(metrics) // len(METHODS),
        "mean_auc_by_method": {method: round(score, 4) for method, score in summary.roc_auc.items()},
    })


if __name__ == "__main__":
    main()
