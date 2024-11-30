import unittest
from cfutils.run import report_mutation

class TestRunFunc(unittest.TestCase):
    """Test cfutils/run.py."""

    def test_report_mutation(self):
        """Test report_mutation function."""
        report_mutation(
            query_ab1_file="./data/B5-M13R_B07.ab1",
            subject_fasta_file="./data/ref.fa",
            output_dir="./temp",
            file_basename="test",
            report_all_sites=True,
            report_mut_plot=False,
        )
        self.assertTrue(True)  # Add more specific assertions as needed

if __name__ == "__main__":
    unittest.main() 