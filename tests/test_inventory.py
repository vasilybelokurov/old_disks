import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
from build_inventory import union  # noqa: E402


def test_union_removes_duplicates_within_a_source():
    """The same star listed twice by one source (e.g. Walker+2023 Hectochelle and M2FS rows) is counted once."""
    m = pd.DataFrame(dict(galaxy="Draco", source="Walker+2023", ra=[260.0, 260.0 + 1e-5, 260.1],
                          dec=[57.9, 57.9, 57.9], feh=[-2.0, -2.1, -1.5]))
    assert union(m).in_union.sum() == 2
