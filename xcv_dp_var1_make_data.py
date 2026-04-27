#!/usr/bin/env python
"""Write data for changepoints::CV.search.DP.VAR1 comparisons."""

from __future__ import annotations

import sys
from pathlib import Path

from xdp_var1_make_data import main as write_base


def main() -> None:
    data_file = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xcv_dp_var1_data.txt")
    sys.argv = [sys.argv[0], str(data_file)]
    write_base()


if __name__ == "__main__":
    main()
