#!/usr/bin/env python3
"""Export gene-specific, precomputed AlphaGenome Atlas RNA scores.

Install: python -m pip install -U alphagenome
Example (100 bp within human SCN1A P1a, GRCh38):
  python download_atlas_rna.py --gene SCN1A --region chr2:166149101-166149200

In a Colab code cell, after installing alphagenome and uploading this file:
  %run download_atlas_rna.py --metadata-only
  %run download_atlas_rna.py --secret-name YOUR_SECRET_NAME --gene SCN1A
The default Colab Secret name is ALPHA_GENOME_API_KEY. Enable Notebook access
for that Secret. Use %run to read Secrets directly through the notebook kernel.
An existing ALPHAGENOME_API_KEY or ALPHA_GENOME_API_KEY environment variable
takes priority; local execution can also use hidden interactive input.

Region input is 1-based inclusive. The SDK receives 0-based half-open bounds.
API access has not been exercised in the authoring session (no API key).
Source contract checked 2026-11-03:
https://github.com/google-deepmind/alphagenome/blob/main/src/alphagenome/atlas/atlas.py
"""

import argparse
import csv
from datetime import datetime, timezone
import getpass
import gzip
from importlib.metadata import version
import json
import os
from pathlib import Path
import re


def load_api_key(secret_name="ALPHA_GENOME_API_KEY"):
    """Read credentials from the environment, Colab Secrets, or a local prompt."""
    for name in ("ALPHAGENOME_API_KEY", "ALPHA_GENOME_API_KEY"):
        value = os.environ.get(name, "").strip()
        if value:
            return value

    try:
        from google.colab import userdata
    except ImportError:
        userdata = None

    if userdata is not None:
        try:
            key = userdata.get(secret_name)
        except userdata.SecretNotFoundError:
            raise SystemExit(
                f"Colab Secret {secret_name!r} was not found. Use --secret-name "
                "with the exact name shown in Colab's Secrets panel."
            ) from None
        except userdata.NotebookAccessError:
            raise SystemExit(
                f"Enable Notebook access for {secret_name!r} in Colab's Secrets panel."
            ) from None
        except Exception:
            raise SystemExit(
                "Could not read Colab Secrets. Run this script with %run in a "
                "Colab code cell. For !python execution, first read the Secret "
                "in a notebook cell and set os.environ['ALPHAGENOME_API_KEY']."
            ) from None
        if not isinstance(key, str) or not key.strip():
            raise SystemExit(f"Colab Secret {secret_name!r} has an empty API key.")
        return key.strip()

    key = getpass.getpass("AlphaGenome API key (hidden input): ").strip()
    if not key:
        raise SystemExit("An API key is required. No key is saved to disk.")
    return key


def parse_region(region):
    match = re.fullmatch(r"(chr[0-9XYM]+):([0-9,]+)-([0-9,]+)", region)
    if not match:
        raise ValueError("Use chr2:166149101-166149200 (1-based inclusive).")
    chrom, first, last = match.groups()
    first, last = int(first.replace(",", "")), int(last.replace(",", ""))
    if first < 1 or last < first:
        raise ValueError("Region requires 1 <= start <= end.")
    return chrom, first - 1, last


def export_scores(adata, output_dir):
    """Keep every returned gene/variant row and the original track ordering."""
    import numpy as np
    import pandas as pd

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs.copy()
    tracks = adata.var.copy()
    for column in ("variant", "interval"):
        if column in obs:
            obs[column] = obs[column].map(str)
    obs.index.name = "score_row_id"
    tracks.index.name = "track_id"
    obs.to_csv(output_dir / "rows.csv")
    tracks.to_csv(output_dir / "tracks.csv")

    def dense(matrix):
        return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)

    raw = dense(adata.X)
    quantile_layer = adata.layers.get("quantiles")
    quantiles = dense(quantile_layer) if quantile_layer is not None else None
    if raw.shape != (len(obs), len(tracks)):
        raise ValueError("Score dimensions do not match row and track metadata.")
    pd.DataFrame(raw, index=obs.index, columns=tracks.index).to_csv(
        output_dir / "raw_scores.csv.gz", compression="gzip", index_label="score_row_id"
    )
    if quantiles is not None:
        pd.DataFrame(quantiles, index=obs.index, columns=tracks.index).to_csv(
            output_dir / "quantile_scores.csv.gz",
            compression="gzip", index_label="score_row_id"
        )

    obs_columns = list(obs.columns)
    track_columns = list(tracks.columns)
    with gzip.open(output_dir / "scores.long.csv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["score_row_id", "track_id"]
            + ["row_" + str(c) for c in obs_columns]
            + ["track_" + str(c) for c in track_columns]
            + ["raw_score", "quantile_score"]
        )
        row_values = list(obs.itertuples(index=False, name=None))
        track_values = list(tracks.itertuples(index=False, name=None))
        for i, row in enumerate(row_values):
            for j, track in enumerate(track_values):
                writer.writerow(
                    [obs.index[i], tracks.index[j]] + list(row) + list(track)
                    + [raw[i, j], quantiles[i, j] if quantiles is not None else ""]
                )
    return {"score_rows": len(obs), "tracks": len(tracks),
            "has_quantiles": quantiles is not None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--secret-name", default="ALPHA_GENOME_API_KEY",
                        help="Colab Secret name (default: ALPHA_GENOME_API_KEY).")
    parser.add_argument("--gene", default="SCN1A", help="Target human gene symbol.")
    parser.add_argument("--region", default="chr2:166149101-166149200")
    parser.add_argument("--ontology", action="append", default=None,
                        help="Exact ontology_curie from metadata; repeat to include more.")
    parser.add_argument("--scorer", choices=["RNA_SEQ", "RNA_SEQ_ACTIVE"],
                        default="RNA_SEQ")
    parser.add_argument("--metadata-only", action="store_true")
    parser.add_argument("--out", type=Path, default=Path("atlas_rna_export"))
    args = parser.parse_args()
    chrom, start0, end0 = parse_region(args.region)
    from alphagenome.atlas import atlas
    from alphagenome.data import genome

    key = load_api_key(args.secret_name)
    client = atlas.create(api_key=key, timeout=60)
    metadata = client.scorer_metadata()
    if args.scorer not in metadata:
        raise SystemExit(f"Scorer {args.scorer} absent. Available: {sorted(metadata)}")
    args.out.mkdir(parents=True, exist_ok=True)
    track_metadata = metadata[args.scorer].track_metadata
    track_metadata.to_csv(args.out / "available_tracks.csv", index_label="track_id")
    print(f"Saved available track metadata: {len(track_metadata)} tracks")
    if args.metadata_only:
        return
    if args.ontology:
        available = set(track_metadata.get("ontology_curie", []))
        missing = set(args.ontology) - available
        if missing:
            raise SystemExit(f"Ontology terms absent from this scorer: {sorted(missing)}")
    result = client.query_interval(
        interval=genome.Interval(chrom, start0, end0),
        requested_scorers=[args.scorer],
        gene_names=[args.gene],
        ontology_terms=args.ontology,
    )
    if args.scorer not in result or result[args.scorer].n_obs == 0:
        raise SystemExit("No scores returned for this gene/region/filter; this is not a zero effect.")
    adata = result[args.scorer]
    if "gene_name" in adata.obs:
        returned_genes = set(adata.obs["gene_name"].dropna())
        if returned_genes != {args.gene}:
            raise RuntimeError(f"Unexpected returned genes: {returned_genes}")
    summary = export_scores(adata, args.out)
    manifest = {
        "retrieved_utc": datetime.now(timezone.utc).isoformat(),
        "alphagenome_package_version": version("alphagenome"),
        "assembly": "GRCh38/hg38", "annotation_in_atlas_report": "GENCODE v46",
        "gene_filter": args.gene, "variant_region_1based_inclusive": args.region,
        "variant_region_0based_halfopen": [chrom, start0, end0],
        "ontology_filter": args.ontology, "scorer": args.scorer,
        "meaning": "Natural-log ALT/REF fold change over gene exon mask" if args.scorer == "RNA_SEQ"
                   else "Maximum REF/ALT mean exon signal; not fold change",
        "query_interval_is_variant_selection_not_rna_aggregation_window": True,
        "values_are_precomputed_predictions_not_experimental_rna_seq_reads": True,
        **summary,
    }
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Saved {summary['score_rows']} gene/variant rows x {summary['tracks']} tracks to {args.out}")


if __name__ == "__main__":
    main()
