#!/usr/bin/env python3

def infer_strandedness(file_path):
    forward_fraction = None
    reverse_fraction = None

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            # Paired-end keys
            if 'Fraction of reads explained by "1++,1--,2+-,2-+"' in line:
                forward_fraction = float(line.split(":")[-1].strip())
            elif 'Fraction of reads explained by "1+-,1-+,2++,2--"' in line:
                reverse_fraction = float(line.split(":")[-1].strip())

            # Single-end keys
            elif 'Fraction of reads explained by "++,--"' in line:
                forward_fraction = float(line.split(":")[-1].strip())
            elif 'Fraction of reads explained by "+-,-+"' in line:
                reverse_fraction = float(line.split(":")[-1].strip())

    # Handle missing values gracefully
    if forward_fraction is None or reverse_fraction is None:
        print("0")   # default to unstranded
        return

    # Determine strandedness
    if forward_fraction > 0.2 and reverse_fraction < 0.2:
        strandedness = "1"  # Forward stranded
    elif reverse_fraction > 0.2 and forward_fraction < 0.2:
        strandedness = "2"  # Reverse stranded
    else:
        strandedness = "0"  # Unstranded / mixed

    print(strandedness)


if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python infer_strandedness.py <rseqc_output_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    infer_strandedness(file_path)
