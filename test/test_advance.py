#!/usr/bin/env python3
"""Unit tests for cfutils advanced mutation reporting."""

import unittest
import tempfile
import os
from cfutils.run import report_mutation

class TestFunc(unittest.TestCase):
    """Test advanced mutation reporting in cfutils."""

    def test_plot_mutation(self) -> None:
        """Test report_mutation with plot output enabled and check output files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                report_mutation(
                    query_ab1_file="./data/B5-M13R_B07.ab1",
                    subject_fasta_file="./data/ref.fa",
                    output_dir=tmpdir,
                    file_basename="test",
                    report_all_sites=True,
                    report_mut_plot=True,
                )
                pdf_path = os.path.join(tmpdir, "test.pdf")
                tsv_path = os.path.join(tmpdir, "test.tsv")
                if not os.path.exists(tsv_path):
                    print(f"Temp dir contents: {os.listdir(tmpdir)}")
                self.assertTrue(os.path.exists(tsv_path), f"Missing TSV: {tsv_path}")
                if os.path.exists(tsv_path):
                    with open(tsv_path) as f:
                        lines = f.readlines()
                    # If there are mutations, PDF should exist
                    if len(lines) > 1:
                        self.assertTrue(os.path.exists(pdf_path), f"Missing PDF: {pdf_path} (but TSV has mutations)")
                    else:
                        if not os.path.exists(pdf_path):
                            print(f"Warning: No mutations found, so PDF was not generated. TSV content: {lines}")
            except Exception as e:
                self.fail(f"report_mutation raised an exception: {e}")

if __name__ == "__main__":
    unittest.main()
