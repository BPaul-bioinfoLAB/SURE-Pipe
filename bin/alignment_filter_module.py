#!/usr/bin/env python3
import sys
import os
import pandas as pd
import argparse
import operator

# Normal operator mapping
ops = {
    ">": operator.gt,
    "<": operator.lt,
    ">=": operator.ge,
    "<=": operator.le,
    "==": operator.eq,
    "!=": operator.ne
}

# Inverse operator mapping (for unique regions)
inverse_ops = {
    ">": operator.lt,
    "<": operator.gt,
    ">=": operator.le,
    "<=": operator.ge,
}
def apply_filter(series, values, op, interval_mode="inside"):
    """
    Flexible filter:
    - If values = [x], applies operator (ops[op])
    - If values = [x, y], applies interval filter
    """
    if not isinstance(values, (list, tuple)):
        values = [values]  # wrap single floats into list

    if len(values) == 2:  # interval mode
        if interval_mode == "inside":
            return series.between(values[0], values[1])
        elif interval_mode == "outside":
            return ~series.between(values[0], values[1])
        else:
            raise ValueError("interval_mode must be 'inside' or 'outside'")
    elif len(values) == 1:  # single threshold mode
        return ops[op](series, values[0])
    else:
        raise ValueError("Invalid number of values passed: must be 1 or 2")

def main():
    parser = argparse.ArgumentParser(description="Parse aligned core-genome file and filter aligned & unaligned core-genome with flexible operators.")

    # Input/output
    parser.add_argument("-i", "--input", required=True, help="Path to the inter-BLAST output file")
    parser.add_argument("-ib", "--ibed", required=True, help="unaligned core-genome BED file")
    parser.add_argument("-o", "--output_dir", default=".", help="Directory to save output files")

    # Unique region thresholds (will use inverse ops internally)
    parser.add_argument("-uqc", "--unique_qcov", type=float, default=0, help="Query coverage % for Unique genomic region (default=0)")
    parser.add_argument("-uqc_op", "--unique_qcov_op", choices=ops.keys(), default=">", help="Operator for qcovs filter in unique regions (will be inverted)")
    parser.add_argument("-uid", "--unique_ident", type=float, default=85.0, help="Identity % for Unique genomic region (default=85.0)")
    parser.add_argument("-uid_op", "--unique_ident_op", choices=ops.keys(), default=">", help="Operator for pident filter in unique regions (will be inverted)")

    # shared region thresholds (use as-is)
    parser.add_argument("-cqc_op", "--shared_qcov_op", choices=ops.keys(), default=">", help="Operator for qcovs filter in shared regions")
    parser.add_argument("-cqc_mode", "--shared_qcov_mode", choices=["inside", "outside"], default="inside", help="Range mode for shared qcov (inside/outside, default=inside)")
    parser.add_argument("-cid_op", "--shared_ident_op", choices=ops.keys(), default="<", help="Operator for pident filter in shared regions")
    parser.add_argument("-cid_mode", "--shared_ident_mode", choices=["inside", "outside"], default="inside", help="Range mode for shared identity (inside/outside, default=inside)")
    
    def parse_range(value):
    # Convert "80" -> [80.0], "80,95" -> [80.0, 95.0]
        parts = value.split(",")
        return [float(p) for p in parts]

    parser.add_argument("-cqc", "--shared_qcov", type=parse_range, required=True, help="Query coverage % for shared genomic region (single or range, e.g. 80 or 80,95)")
    parser.add_argument("-cid", "--shared_ident", type=parse_range, required=True, help="Identity % for shared genomic region (single or range, e.g. 90 or 90,100)")

    args = parser.parse_args()
    
    #Step 1:Ensure output dir exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Define output file names
    output_valid = os.path.join(args.output_dir, "core_vs_neighbour_shared_regions.csv")
    output_trimmed = os.path.join(args.output_dir, "core_vs_neighbour_unique_regions.tsv")
    output_bed = os.path.join(args.output_dir, "core_unaligned_coordinates.bed")
    open(output_bed, 'w').close()

    print(f"\n🔄 Loading file: {args.input}...")

    #Step 2:Load BLAST output
    columns = ["qseqid", "sseqid", "stitle", "qcovs", "pident", "qstart", "qend", "sstart", "send", 
               "length", "evalue", "bitscore"]
    blast_df = pd.read_csv(args.input, sep="\t", header=None, names=columns)

    # ---------- Unique regions filter (INVERSE ops applied) ----------
    final_df = blast_df[
        inverse_ops[args.unique_qcov_op](blast_df["qcovs"], args.unique_qcov) &
        inverse_ops[args.unique_ident_op](blast_df["pident"], args.unique_ident)
    ]
    final_df.to_csv(output_trimmed, sep="\t", index=False, header=False)
    print(f"✅ Unique (inverted op) BLAST output saved to: {output_trimmed}")

    # ---------- shared regions filter (normal ops) ----------
    shared_df = blast_df[
        apply_filter(blast_df["qcovs"], args.shared_qcov, args.shared_qcov_op, args.shared_qcov_mode) &
        apply_filter(blast_df["pident"], args.shared_ident, args.shared_ident_op, args.shared_ident_mode)
    ]

    shared_filt_df = shared_df[[
        "qseqid", "sseqid", "qcovs", "pident", "qstart", "qend", "sstart", "send", 
        "length", "evalue", "bitscore"
    ]]
    shared_filt_df.to_csv(output_valid, sep="\t", index=False, header=False)
    print(f"✅ shared BLAST output saved to: {output_valid}")

    # ---------- Handle BED file ----------
    if os.path.isfile(args.ibed) and os.path.getsize(args.ibed) > 0:
        bed_df = pd.read_csv(args.ibed, sep="\t", header=None, names=["qseqid", "qstart", "qend", "name"])
        bed_df = bed_df.sort_values(by=["qseqid", "qstart"])
        bed_df.to_csv(output_bed, sep="\t", index=False, header=False)
        print(f"✅ BED file sorted and saved: {output_bed}")
    else:
        print("⚠️ Provided BED file is empty or missing. Skipping BED file creation.")

if __name__ == "__main__":
    main()
