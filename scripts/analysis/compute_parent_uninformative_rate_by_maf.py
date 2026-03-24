#!/usr/bin/env python3

import argparse
import math
import os
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import polars as pl


CHM13_PAR1_REGION = "chrX:1-2394410"
CHM13_PAR2_REGION = "chrX:153925834-154259566"


ROUNDED_MAF_BIN_EDGES: List[float] = [
    0.0,
    0.00021,
    0.00042,
    0.00064,
    0.001,
    0.0016,
    0.0022,
    0.003,
    0.004,
    0.0054,
    0.0072,
    0.0094,
    0.0126,
    0.0172,
    0.0244,
    0.0369,
    0.0601,
    0.1018,
    0.1661,
    0.2556,
    0.3724,
    0.5,
]

ROUNDED_MAF_LABELS: List[str] = [
    "singleton",
    "0.00021-0.00042",
    "0.00042-0.00064",
    "0.00064-0.001",
    "0.001-0.0016",
    "0.0016-0.0022",
    "0.0022-0.003",
    "0.003-0.004",
    "0.004-0.0054",
    "0.0054-0.0072",
    "0.0072-0.0094",
    "0.0094-0.0126",
    "0.0126-0.0172",
    "0.0172-0.0244",
    "0.0244-0.037",
    "0.037-0.06",
    "0.06-0.102",
    "0.102-0.166",
    "0.166-0.255",
    "0.255-0.372",
    "0.372-0.5",
]


@dataclass(frozen=True)
class Trio:
    child: str
    father: str
    mother: str


def _parse_ped_trios(ped_path: str) -> List[Trio]:
    """Parse pedigree-like files and return (child, father, mother) trios.

    This repo uses multiple lightweight formats:
    - 3 columns: child father mother  (e.g., resources/pedigrees/1kgp.ped)
    - 5+ columns: family child father mother sex [phenotype...] (e.g., resources/pedigrees/trios_only.ped)
    - Standard .ped is also accepted: family child father mother sex phenotype [...]
    """

    trios: List[Trio] = []
    with open(ped_path, "r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            # Support comma-separated trio lists (child,father,mother)
            if "," in line:
                parts = [p.strip() for p in line.split(",") if p.strip()]
            else:
                parts = line.split()

            if len(parts) == 3:
                child, father, mother = parts
            elif len(parts) >= 4:
                # family_id, child, father, mother, ...
                child = parts[1]
                father = parts[2]
                mother = parts[3]
            else:
                continue
            if father in ("0", "-9", "NA", ".") or mother in ("0", "-9", "NA", "."):
                continue
            if father == child or mother == child or father == mother:
                continue
            trios.append(Trio(child=child, father=father, mother=mother))
    return trios


def _rounded_maf_label(
    maf: Optional[float],
    minor_allele_count: Optional[int],
) -> Optional[str]:
    if maf is None or math.isnan(maf):
        return None

    if minor_allele_count == 1:
        return "singleton"

    # Match pipeline behavior: bins are (start, end] with start=0 excluded.
    # Here we assume maf is in [0, 0.5].
    for start, end, label in zip(
        ROUNDED_MAF_BIN_EDGES[1:-1],
        ROUNDED_MAF_BIN_EDGES[2:],
        ROUNDED_MAF_LABELS[1:],
    ):
        if maf > start and maf <= end:
            return label

    # If maf is <= first edge and not singleton (e.g. MAC==0 or no MAC), treat as None.
    return None


def _try_get_info_float(info: Dict, key: str) -> Optional[float]:
    if info is None:
        return None
    value = info.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return None
        value = value[0]
    try:
        return float(value)
    except Exception:
        return None


def _try_get_info_int(info: Dict, key: str) -> Optional[int]:
    if info is None:
        return None
    value = info.get(key)
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if len(value) == 0:
            return None
        value = value[0]
    try:
        return int(value)
    except Exception:
        return None


def _get_maf_and_mac(
    info: Dict,
    gt: Optional[np.ndarray],
    prefer_info_maf: bool,
) -> Tuple[Optional[float], Optional[int]]:
    """Return (maf, minor_allele_count).

    Priority:
    1) INFO/MAF (used elsewhere in this repo) if prefer_info_maf.
    2) INFO/AC + INFO/AN -> compute maf and MAC.
    3) INFO/AF + INFO/AN -> approximate AC via round(AF*AN) then compute.
    4) Fallback: compute AC/AN from genotypes (gt) across all samples.
    """

    if prefer_info_maf:
        maf_val = _try_get_info_float(info, "MAF")
        if maf_val is not None:
            maf_val = min(maf_val, 1.0 - maf_val)
            mac = _try_get_info_int(info, "MAC")
            if mac is not None:
                return maf_val, mac
            ac = _try_get_info_int(info, "AC")
            an = _try_get_info_int(info, "AN")
            if ac is not None and an is not None and an > 0:
                return maf_val, min(ac, an - ac)
            return maf_val, None

    ac = _try_get_info_int(info, "AC")
    an = _try_get_info_int(info, "AN")
    if ac is not None and an is not None and an > 0:
        mac = min(ac, an - ac)
        maf = mac / an
        return maf, mac

    af = _try_get_info_float(info, "AF")
    if af is not None and an is not None and an > 0:
        ac_est = int(round(af * an))
        mac = min(ac_est, an - ac_est)
        maf = mac / an
        return maf, mac

    if gt is None:
        return None, None

    # Fallback: compute AC/AN from GT across all samples.
    # Assume biallelic and diploid (PAR1 behaves diploid in both sexes).
    g0 = gt[:, 0]
    g1 = gt[:, 1]
    called = (g0 >= 0) & (g1 >= 0)
    if called.sum() == 0:
        return None, None
    an_calc = int(called.sum() * 2)
    ac_calc = int((g0[called] + g1[called]).sum())
    mac = min(ac_calc, an_calc - ac_calc)
    if an_calc == 0:
        return None, None
    maf = mac / an_calc
    return maf, mac

from cyvcf2 import VCF
from tqdm import tqdm

def extract_gt_types_dense(bcf_path: str) -> tuple[np.ndarray, list[str]]:
    vcf = VCF(bcf_path)
    samples = vcf.samples
    n_samples = len(samples)

    # Preallocate using index-derived count
    n_variants = vcf.num_records

    # 0=homref, 1=het, 2=homalt, 3=missing/unknown
    gt = np.empty((n_variants, n_samples), dtype=np.uint8)

    i = 0
    for rec in tqdm(vcf, total=n_variants):
        # rec.gt_types is already an int array; avoid astype/copy if possible
        gt[i, :] = rec.gt_types
        i += 1

    # Safety in case count metadata was stale (rare)
    if i != n_variants:
        gt = gt[:i, :]

    return gt

def _get_gt_matrix_from_cyvcf2(variant) -> Optional[np.ndarray]:
    """Return numpy array GT of shape (n_samples, 2) with -1 for missing.

    cyvcf2 exposes per-sample genotypes as a list of [a1, a2, phased] rows.
    """

    g = variant.genotypes
    if g is None:
        return None
    gt = np.asarray(g)
    if gt.ndim != 2 or gt.shape[1] < 2:
        return None
    return gt[:, :2]


def _boolean_missing(gt2: np.ndarray) -> np.ndarray:
    return gt2==3


def _boolean_het(gt2: np.ndarray) -> np.ndarray:
    return gt2==1
    return (gt2[:, 0] >= 0) & (gt2[:, 1] >= 0) & (gt2[:, 0] != gt2[:, 1])


def _boolean_hom(gt: np.ndarray) -> np.ndarray:
    """True for called homozygous (0/0 or 1/1) genotypes."""
    ref=gt==0
    alt=gt==2
    return ref | alt
    return (gt2[:, 0] >= 0) & (gt2[:, 1] >= 0) & (gt2[:, 0] == gt2[:, 1])


def compute_uninformative_by_rounded_maf(
    panel_bcf: str,
    ped_path: str,
    region: str,
    genome: str,
    output_parquet: str,
    prefer_info_maf: bool,
    max_variants: Optional[int] = None,
    write_parquet: bool = True,
) -> pl.DataFrame:
    from cyvcf2 import VCF  # local import to keep import errors obvious

    trios_all = _parse_ped_trios(ped_path)
    if len(trios_all) == 0:
        raise ValueError(f"No trios parsed from ped: {ped_path}")

    vcf = VCF(panel_bcf)
    samples = list(vcf.samples)
    sample_to_index = {s: i for i, s in enumerate(samples)}

    trios: List[Trio] = [
        t
        for t in trios_all
        if t.child in sample_to_index and t.father in sample_to_index and t.mother in sample_to_index
    ]

    if len(trios) == 0:
        raise ValueError(
            "No complete trios found in VCF header. "
            f"Parsed {len(trios_all)} from ped, VCF has {len(samples)} samples."
        )

    child_idx = np.fromiter((sample_to_index[t.child] for t in trios), dtype=np.int32)
    dad_idx = np.fromiter((sample_to_index[t.father] for t in trios), dtype=np.int32)
    mom_idx = np.fromiter((sample_to_index[t.mother] for t in trios), dtype=np.int32)

    counts: Dict[str, Dict[str, int]] = {}
    variants_seen = 0

    # Stream region
    gts=extract_gt_types_dense(bcfpath)
    child_gts=gts[:,child_idx]
    dad_gts=gts[:, dad_idx]
    mom_gts=gts[:,]
    for variant in vcf(region):
        # variants_seen += 1
        # if max_variants is not None and variants_seen > max_variants:
        #     break

        gt = _get_gt_matrix_from_cyvcf2(variant)
        if gt is None:
            continue

        info = getattr(variant, "INFO", {}) or {}
        maf, mac = _get_maf_and_mac(info, gt, prefer_info_maf=prefer_info_maf)
        label = _rounded_maf_label(maf=maf, minor_allele_count=mac)
        if label is None:
            continue

        child_gt = gt[child_idx]
        dad_gt = gt[dad_idx]
        mom_gt = gt[mom_idx]

        child_het = _boolean_het(child_gts)

        dad_missing = _boolean_missing(dad_gts)
        mom_missing = _boolean_missing(mom_gts)

        dad_het = _boolean_het(dad_gts)
        mom_het = _boolean_het(mom_gts)

        dad_hom = _boolean_hom(dad_gts)
        mom_hom = _boolean_hom(mom_gts)

        dad_uninf = dad_missing | dad_het
        mom_uninf = mom_missing | mom_het

        # Define uninformative categories among proband heterozygous sites:
        # 1) both parents het
        # 2) both parents missing
        # 3) one parent missing AND the other het
        # 4) informative (at least one parent homozygous-called)
        both_parents_het = child_het & dad_het & mom_het
        both_parents_missing = child_het & dad_missing & mom_missing
        one_missing_one_het = child_het & ((dad_missing & mom_het) | (dad_het & mom_missing))

        # Matches the original long-term definition of "uninformative": both parents are missing or het
        uninformative = child_het & dad_uninf & mom_uninf

        # For diploid biallelic GTs, this is equivalent to "at least one parent is called-homozygous".
        # Using the complement ensures the 4 categories partition proband hets even if unexpected GT
        # encodings appear.
        informative = child_het & ~uninformative

        n_child_het = int(child_het.sum())
        if n_child_het == 0:
            continue

        n_both_het = int(both_parents_het.sum())
        n_both_missing = int(both_parents_missing.sum())
        n_one_missing_one_het = int(one_missing_one_het.sum())
        # By construction, uninformative is the union of the three uninformative subtypes.
        n_uninf = n_both_het + n_both_missing + n_one_missing_one_het
        n_informative = int(informative.sum())

        # Sanity: categories should partition all proband hets
        # (informative is the complement, so this should hold barring overflow/bugs)
        if (n_uninf + n_informative) != n_child_het:
            raise RuntimeError(
                "Category partition failed: "
                f"n_child_het={n_child_het} n_uninf={n_uninf} n_informative={n_informative}"
            )

        if label not in counts:
            counts[label] = {
                "n_proband_hets": 0,
                "n_parent_uninformative": 0,
                "n_uninf_both_parents_het": 0,
                "n_uninf_both_parents_missing": 0,
                "n_uninf_one_missing_one_het": 0,
                "n_informative": 0,
            }

        counts[label]["n_proband_hets"] += n_child_het
        counts[label]["n_parent_uninformative"] += n_uninf
        counts[label]["n_uninf_both_parents_het"] += n_both_het
        counts[label]["n_uninf_both_parents_missing"] += n_both_missing
        counts[label]["n_uninf_one_missing_one_het"] += n_one_missing_one_het
        counts[label]["n_informative"] += n_informative

    rows = []
    for label in ROUNDED_MAF_LABELS:
        if label not in counts:
            continue
        n_hets = counts[label]["n_proband_hets"]
        n_uninf = counts[label]["n_parent_uninformative"]
        n_both_het = counts[label]["n_uninf_both_parents_het"]
        n_both_missing = counts[label]["n_uninf_both_parents_missing"]
        n_one_missing_one_het = counts[label]["n_uninf_one_missing_one_het"]
        n_informative = counts[label]["n_informative"]

        pct = (n_uninf / n_hets * 100.0) if n_hets > 0 else float("nan")
        rows.append(
            {
                "genome": genome,
                "region": region,
                "rounded_MAF": label,
                "n_proband_hets": n_hets,
                "n_parent_uninformative": n_uninf,
                "n_uninf_both_parents_het": n_both_het,
                "n_uninf_both_parents_missing": n_both_missing,
                "n_uninf_one_missing_one_het": n_one_missing_one_het,
                "n_informative": n_informative,
                "pct_parent_uninformative": pct,
                "pct_uninf_both_parents_het": (n_both_het / n_hets * 100.0) if n_hets > 0 else float("nan"),
                "pct_uninf_both_parents_missing": (n_both_missing / n_hets * 100.0) if n_hets > 0 else float("nan"),
                "pct_uninf_one_missing_one_het": (n_one_missing_one_het / n_hets * 100.0) if n_hets > 0 else float("nan"),
                "pct_informative": (n_informative / n_hets * 100.0) if n_hets > 0 else float("nan"),
                "n_trios": len(trios),
                "n_variants_seen": variants_seen,
                "panel_bcf": os.path.basename(panel_bcf),
            }
        )

    df = pl.DataFrame(rows)

    if df.height == 0:
        raise ValueError(
            "No rows produced. This likely means variants lacked INFO/MAF/AC/AN, "
            "or the region has no variants, or GT extraction failed."
        )

    if write_parquet:
        df.write_parquet(output_parquet)
    return df


def _bcf_path_from_template(template: str, replacement: str, token: str = "PAR1") -> str:
    if token not in template:
        raise ValueError(
            f"panel BCF template does not contain token '{token}': {template}. "
            "Provide a template path that includes the region token so it can be replaced."
        )
    return template.replace(token, replacement)


def compute_uninformative_chr1_22_and_par2_from_template(
    panel_bcf_template: str,
    ped_path: str,
    genome: str,
    output_parquet: str,
    prefer_info_maf: bool,
    max_variants: Optional[int] = None,
    template_token: str = "PAR1",
) -> pl.DataFrame:
    """Compute stats for chr1-chr22 and PAR2, using a BCF path template.

    The repo's testing runs often use per-region working directories/BCFs where
    the only difference is the region token (e.g. PAR1 -> chr18, PAR1 -> PAR2).
    This helper derives each region-specific BCF path by string replacement.

    Missing BCFs are skipped (with a warning to stderr).
    """

    dfs: List[pl.DataFrame] = []
    missing: List[str] = []

    for chrom in range(1, 23):
        chrom_label = f"chr{chrom}"
        panel_bcf = _bcf_path_from_template(
            panel_bcf_template, chrom_label, token=template_token
        )
        if not os.path.exists(panel_bcf):
            missing.append(panel_bcf)
            continue
        dfs.append(
            compute_uninformative_by_rounded_maf(
                panel_bcf=panel_bcf,
                ped_path=ped_path,
                region=chrom_label,
                genome=genome,
                output_parquet=output_parquet,
                prefer_info_maf=prefer_info_maf,
                max_variants=max_variants,
                write_parquet=False,
            )
        )

    # PAR2
    panel_bcf_par2 = _bcf_path_from_template(
        panel_bcf_template, "PAR2", token=template_token
    )
    if os.path.exists(panel_bcf_par2):
        dfs.append(
            compute_uninformative_by_rounded_maf(
                panel_bcf=panel_bcf_par2,
                ped_path=ped_path,
                region=CHM13_PAR2_REGION,
                genome=genome,
                output_parquet=output_parquet,
                prefer_info_maf=prefer_info_maf,
                max_variants=max_variants,
                write_parquet=False,
            )
        )
    else:
        missing.append(panel_bcf_par2)

    if not dfs:
        raise ValueError(
            "No region BCFs were found from the provided template. "
            "Check the template path and token replacement logic."
        )

    if missing:
        print(
            f"Warning: {len(missing)} region BCF(s) missing; skipped. First few:\n"
            + "\n".join(missing[:10]),
            file=sys.stderr,
        )

    out = pl.concat(dfs, how="vertical")
    out.write_parquet(output_parquet)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute % of proband heterozygous sites where both parents are uninformative "
            "(each parent missing OR heterozygous), stratified by existing rounded_MAF bins. "
            "Intended for CHM13v2.0; start with PAR1."
        )
    )

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--panel-bcf",
        help=(
            "Single-region biallelic BCF/VCF with trio samples. "
            "Example: .../1KGP.CHM13v2.0.PAR1.snp_indel.phasing_qual_pass.biallelic.bcf"
        ),
    )
    input_group.add_argument(
        "--panel-bcf-template",
        help=(
            "Template path used for multi-region mode. The script will derive region-specific BCFs "
            "by replacing the token (default 'PAR1') with 'chr1'..'chr22' and 'PAR2'. "
            "Example template: /.../PAR1_working_.../1KGP.CHM13v2.0.PAR1.snp_indel.phasing_qual_pass.biallelic.bcf"
        ),
    )
    parser.add_argument(
        "--ped",
        default="resources/pedigrees/trios_only.ped",
        help="Pedigree/trios file (default: resources/pedigrees/trios_only.ped)",
    )
    parser.add_argument(
        "--region",
        default=CHM13_PAR1_REGION,
        help=f"Region string (default: {CHM13_PAR1_REGION})",
    )
    parser.add_argument(
        "--chr1-22-and-par2",
        action="store_true",
        help=(
            "Run multi-region mode (chr1-chr22 + PAR2). Requires --panel-bcf-template. "
            "BCF paths are derived by replacing the template token (default 'PAR1')."
        ),
    )
    parser.add_argument(
        "--template-token",
        default="PAR1",
        help="Token to replace in --panel-bcf-template (default: PAR1)",
    )
    parser.add_argument(
        "--genome",
        default="CHM13v2.0",
        help="Genome label to record in output (default: CHM13v2.0)",
    )
    parser.add_argument(
        "--out-parquet",
        default="intermediate_data/proband_parent_uninformative_by_MAF.parquet",
        help="Output parquet path (default: intermediate_data/proband_parent_uninformative_by_MAF.parquet)",
    )
    parser.add_argument(
        "--prefer-info-maf",
        action="store_true",
        help="Prefer INFO/MAF when present (recommended for repo-produced panels).",
    )
    parser.add_argument(
        "--max-variants",
        type=int,
        default=None,
        help="Optional cap on number of variants processed (debug).",
    )

    args = parser.parse_args()

    df = compute_uninformative_by_rounded_maf(
        panel_bcf=args.panel_bcf,
        ped_path=args.ped,
        region=args.region,
        genome=args.genome,
        output_parquet=args.out_parquet,
        prefer_info_maf=args.prefer_info_maf,
        max_variants=args.max_variants,
    )

    # Minimal stdout summary for sanity checks
    with pl.Config(tbl_rows=40, tbl_cols=12):
        print(df.sort("rounded_MAF"))


if __name__ == "__main__":
    main()
