#!/usr/bin/env python
from __future__ import annotations

import sys
from pathlib import Path


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("xmcp_demo_data.txt")
    out.write_text("# populated by xmcp_demo_file.R from mcp::demo_fit\n", encoding="ascii")
    print(f"prepared placeholder {out.name}")
    print("reference source = mcp::demo_fit package data")


if __name__ == "__main__":
    main()
