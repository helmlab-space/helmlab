"""CIEDE2000 of hl.metric.ciede2000 against the 34 pairs of Sharma, Wu & Dalal (2005).

Color Research & Application 30(1), 21-30. The file tests/data/sharma2005_ciede2000testdata.txt
is ciede2000testdata.txt from https://hajim.rochester.edu/ece/sites/gsharma/ciede2000/,
byte for byte (L1, a1, b1, L2, a2, b2, dE00 to 4 decimals). Each value rounded to 4 decimals
must equal the printed one, in both orders of the two colours.
"""
import hashlib
from pathlib import Path

import numpy as np

from helmlab.helmlab import _ciede2000

DATA = Path(__file__).parent / "data" / "sharma2005_ciede2000testdata.txt"
SHA256 = "44aebb39107128328add54fbef5ac8ee89909e50508f448a1580adea2058a4b8"


def test_sharma_2005_all_34_pairs():
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == SHA256
    rows = np.loadtxt(DATA, delimiter="\t")
    assert rows.shape == (34, 7)
    for i, r in enumerate(rows, start=1):
        lab1, lab2, printed = r[0:3], r[3:6], r[6]
        assert round(_ciede2000(lab1, lab2), 4) == printed, f"pair {i}"
        assert round(_ciede2000(lab2, lab1), 4) == printed, f"pair {i}, colours swapped"
