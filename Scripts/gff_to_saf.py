#!/usr/bin/env python3

import argparse


def parse_attributes(attr):
    """
    Parse GFF3 attributes exactly enough for ID extraction.
    No duplicate handling.
    No artificial ID changes.
    """
    d = {}
    for part in attr.strip().split(";"):
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
    return d


def get_id(attrs, id_source):
    """
    Extract selected ID exactly from GFF attribute.

    Examples:
    --id-source Parent  -> gene-slr0612
    --id-source ID      -> cds-BAA16580.1
    --id-source gene    -> slr0612
    --id-source Name    -> BAA16580.1
    --id-source locus_tag -> slr0612
    """
    if id_source == "Parent":
        parent = attrs.get("Parent")
        if not parent:
            return None
        return parent.split(",")[0]

    return attrs.get(id_source)


def main():
    parser = argparse.ArgumentParser(
        description="Convert GFF3 to SAF for featureCounts"
    )

    parser.add_argument("-i", "--gff", required=True, help="Input GFF3 file")
    parser.add_argument("-o", "--saf", required=True, help="Output SAF file")

    parser.add_argument(
        "--feature",
        default="CDS",
        help="Feature type to extract. Default: CDS"
    )

    parser.add_argument(
        "--id-source",
        default="Parent",
        help=(
            "GFF attribute to use as GeneID. "
            "Examples: Parent, ID, gene, Name, locus_tag, protein_id. "
            "Default: Parent"
        )
    )

    args = parser.parse_args()

    with open(args.gff) as fin, open(args.saf, "w") as out:
        out.write("GeneID\tChr\tStart\tEnd\tStrand\n")

        for line in fin:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue

            chrom = fields[0]
            feature_type = fields[2]
            start = fields[3]
            end = fields[4]
            strand = fields[6]
            attr = fields[8]

            if feature_type != args.feature:
                continue

            if strand not in {"+", "-"}:
                continue

            attrs = parse_attributes(attr)
            gene_id = get_id(attrs, args.id_source)

            if not gene_id:
                continue

            out.write(f"{gene_id}\t{chrom}\t{start}\t{end}\t{strand}\n")


if __name__ == "__main__":
    main()
