"""
Stress-test Python and Fortran ArgpCpd implementations on multiple simulated files.
"""

from __future__ import annotations

import ast
import subprocess
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import changepoint as cpt
import numpy as np

from xargpcpd_make_data import simulate_argpcpd_signal


@dataclass(frozen=True)
class StressCase:
    name: str
    seed: int
    n: int
    true_bkps: list[int]
    means: list[float]
    sds: list[float]


CASES = [
    StressCase("base", 31, 600, [0, 200, 400], [0.0, 4.0, -3.0], [0.4, 0.4, 0.5]),
    StressCase("shifted", 41, 600, [0, 180, 390], [1.0, -3.5, 2.5], [0.5, 0.45, 0.55]),
    StressCase("higher_var", 77, 600, [0, 210, 420], [0.0, 5.0, -2.0], [0.6, 0.7, 0.8]),
    StressCase("closer_regimes", 103, 600, [0, 220, 410], [0.0, 2.8, -2.5], [0.35, 0.35, 0.45]),
    StressCase("late_change", 151, 600, [0, 240, 470], [0.5, 4.5, -3.2], [0.45, 0.45, 0.55]),
]


def fit_python(signal: np.ndarray) -> list[int]:
    argp = cpt.ArgpCpd(
        logistic_hazard_h=-5.0,
        scale=3.0,
        noise_level=0.01,
        max_lag=12,
    )
    rs = [argp.step(float(x)) for x in signal]
    return list(cpt.map_changepoints(rs))


def write_signal(path: Path, case: StressCase, signal: np.ndarray) -> None:
    with path.open("w", encoding="ascii") as f:
        f.write(f"# case = {case.name}\n")
        f.write(f"# seed = {case.seed}\n")
        f.write(f"# n = {case.n}\n")
        f.write("# true_bkps = " + " ".join(str(x) for x in case.true_bkps) + "\n")
        f.write("# means = " + " ".join(f"{x:.6f}" for x in case.means) + "\n")
        f.write("# sds = " + " ".join(f"{x:.6f}" for x in case.sds) + "\n")
        np.savetxt(f, signal, fmt="%.6f")


def fit_fortran(data_file: Path) -> list[int]:
    proc = subprocess.run(
        ["xsim_argpcpd_gp_file.exe", str(data_file), "rust"],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in proc.stdout.splitlines():
        if "estimated" in line:
            return [int(tok) for tok in line.split("=")[1].split()]
    raise RuntimeError("Fortran output did not contain an estimated line")


def main() -> None:
    t0 = time.perf_counter()
    exe = Path("xsim_argpcpd_gp_file.exe")
    if not exe.exists():
        raise SystemExit("Build xsim_argpcpd_gp_file.exe first with: make -f Makefile.xcorr xsim_argpcpd_gp_file")

    n_match = 0
    with tempfile.TemporaryDirectory(prefix="xargpcpd_stress_", dir=".") as tmpdir:
        tmpdir_path = Path(tmpdir)
        for case in CASES:
            signal = simulate_argpcpd_signal(case.n, case.true_bkps, case.means, case.sds, case.seed)
            data_file = tmpdir_path / f"{case.name}.txt"
            write_signal(data_file, case, signal)
            py_cps = fit_python(signal)
            ft_cps = fit_fortran(data_file)
            matched = py_cps == ft_cps
            if matched:
                n_match += 1
            print(f"{case.name:14s} match={str(matched):5s} python={py_cps} fortran={ft_cps}")

    print(f"matches = {n_match}/{len(CASES)}")
    print(f"elapsed seconds = {time.perf_counter() - t0:.3f}")


if __name__ == "__main__":
    main()
