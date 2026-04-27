#!/usr/bin/env python

from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_sigma_data.txt")
    out.write_text("0 0\n", encoding="ascii")
    print(f"prepared placeholder {out}")
    print("reference source = repo mcp variance example")


if __name__ == "__main__":
    main()
