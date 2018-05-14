#!/usr/bin/env python3

"""test for cfutils"""

import sys

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError('cfutils requires Python 3.6+!')
