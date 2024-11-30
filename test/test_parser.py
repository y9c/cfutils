import unittest
from cfutils.parser import parse_abi, parse_fasta

class TestParserFunc(unittest.TestCase):
    """Test cfutils/parser.py."""

    def test_parse_abi(self):
        """Test parse_abi function."""
        record = parse_abi("./data/B5-M13R_B07.ab1")
        self.assertIsNotNone(record)

    def test_parse_fasta(self):
        """Test parse_fasta function."""
        record = parse_fasta("./data/ref.fa")
        self.assertIsNotNone(record)

if __name__ == "__main__":
    unittest.main() 