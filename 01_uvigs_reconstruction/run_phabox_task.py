#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run a single PhaBOX2 task with explicit proteins.")
    parser.add_argument("--task", required=True, choices=["phagcn", "phatyp"])
    parser.add_argument("--dbdir", required=True)
    parser.add_argument("--outpth", required=True)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--proteins", required=True)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--len", dest="len", type=int, default=2000)
    parser.add_argument("--midfolder", default="midfolder")
    parser.add_argument("--aai", type=float, default=75.0)
    parser.add_argument("--share", type=float, default=15.0)
    parser.add_argument("--pcov", type=float, default=80.0)
    parser.add_argument("--draw", default="N")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    env_bin = str(Path(sys.executable).resolve().parent)
    os.environ["PATH"] = f"{env_bin}:{os.environ.get('PATH', '')}"

    from phabox2 import phagcn, phatyp

    Path(args.outpth).mkdir(parents=True, exist_ok=True)
    if args.task == "phagcn":
        phagcn.run(args)
    elif args.task == "phatyp":
        phatyp.run(args)
    else:
        raise ValueError(f"Unsupported task: {args.task}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
