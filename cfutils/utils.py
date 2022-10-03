#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-05-27 21:19

"""shared function for package."""

import logging
import sys

try:
    assert sys.version_info > (3, 6)
except AssertionError:
    raise RuntimeError("cfutils requires Python 3.6+!")


def get_logger(name: str) -> logging.Logger:
    """global logging."""
    logger: logging.Logger = logging.getLogger(name)
    if not logger.handlers:
        handler: logging.StreamHandler = logging.StreamHandler()
        formatter: logging.Formatter = logging.Formatter(
            "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        #  logger.setLevel(logging.DEBUG)
        logger.setLevel(logging.INFO)
    return logger


LOGGER: logging.Logger = get_logger(__name__)


def evenchunks(string, chunksize=10):
    out = []
    for i in range(0, len(string), chunksize):
        end = i + chunksize
        out.append(string[i:end])
    return out


def chunked_lines(string, chunksize=10, chunks_per_line=5, spacer=" "):
    chunks = evenchunks(string, chunksize)
    lines = []
    while chunks:
        lines.append(spacer.join(chunks[:chunks_per_line]))
        del chunks[:chunks_per_line]
    return lines
