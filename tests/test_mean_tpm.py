import tempfile
import unittest
from pathlib import Path

import pandas as pd

from Scripts.mean_tpm_calc import mean_tpm_calc


class MeanTpmTests(unittest.TestCase):
    def test_sample_and_condition_underscores_are_unambiguous(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            samplesheet = tmp_path / "samples.csv"
            samplesheet.write_text(
                "sample_id,condition,fastq_1,fastq_2\n"
                "sample_one,control_group,r1.fastq,\n"
                "sample_two,control_group,r2.fastq,\n"
            )

            base = {
                "Geneid": ["gene_a", "gene_b"],
                "Chr": ["chr", "chr"],
                "Start": [1, 11],
                "End": [10, 20],
                "Strand": ["+", "-"],
                "Length": [10, 10],
            }
            files = []
            for sample_id, values in (
                ("sample_one", [10.0, 20.0]),
                ("sample_two", [30.0, 40.0]),
            ):
                frame = pd.DataFrame(base)
                frame[f"{sample_id}_control_group_TPM"] = values
                path = tmp_path / f"{sample_id}_normalized.results"
                frame.to_csv(path, index=False)
                files.append(path)

            output = tmp_path / "Expression.csv"
            mean_tpm_calc(files, samplesheet, output)

            result = pd.read_csv(output)
            self.assertEqual(result["control_group"].tolist(), [20.0, 30.0])


if __name__ == "__main__":
    unittest.main()
