#!/usr/bin/env python3
"""test for cfutils"""

import sys
import unittest


class TestBasicFunc(unittest.TestCase):
    """Test basic"""

    def test_basic(self):
        """Test basic"""
        self.assertGreaterEqual(sys.version_info, (3, 6))


if __name__ == "__main__":
    unittest.main()
