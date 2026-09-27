"""Load precomputed folded-instrument arithmetic for the offline study."""

from __future__ import annotations

import ctypes
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from benchmark_wan_scaling import compact_plan
from sahay_gadget import Bound, SahayGadget
from study_wan_scaling import REFERENCE, paired_masks


class SahayKernel:
    def __init__(self, gadget: SahayGadget, reference: Path, directory: Path):
        directory.mkdir(parents=True, exist_ok=False)
        rank, leaves = gadget.rank, len(gadget.masks)
        (directory / "metadata.txt").write_text(f"{rank}\n{leaves}\n")
        (directory / "amplitude.txt").write_text(compact_plan(gadget.masks, rank))
        for k in range(rank + 1):
            (directory / f"marginal_{k}.txt").write_text(
                compact_plan(paired_masks(gadget.masks, rank, k), rank + k)
            )
        relative = "tools/profile/gadget_contraction_kernel.h"
        header = reference / relative
        if header.read_bytes() != subprocess.check_output(
            ["git", "show", REFERENCE + ":" + relative]
        ):
            raise ValueError("Unpinned contraction kernel")
        source = Path(__file__).with_name("sahay_contraction_kernel.cc")
        library = directory / "kernel.so"
        subprocess.run(
            [
                "c++",
                "-std=c++20",
                "-O3",
                "-march=native",
                "-fno-fast-math",
                "-shared",
                "-fPIC",
                "-I" + str(header.parent),
                "-I" + str(source.resolve().parents[2] / "src"),
                str(source),
                "-o",
                str(library),
            ],
            check=True,
        )
        self.library = ctypes.CDLL(str(library))
        array = np.ctypeslib.ndpointer(dtype=np.complex128, flags="C_CONTIGUOUS")
        bits = np.ctypeslib.ndpointer(dtype=np.uint8, flags="C_CONTIGUOUS")
        self.library.sahay_create.argtypes = [ctypes.c_char_p]
        self.library.sahay_create.restype = ctypes.c_void_p
        self.library.sahay_destroy.argtypes = [ctypes.c_void_p]
        self.library.sahay_marginal.argtypes = [
            ctypes.c_void_p,
            array,
            array,
            bits,
            ctypes.c_size_t,
        ]
        self.library.sahay_marginal.restype = ctypes.c_double
        self.library.sahay_amplitudes.argtypes = [ctypes.c_void_p, array, array, bits, array]
        self.context = self.library.sahay_create(str(directory).encode())
        if not self.context:
            raise ValueError("Invalid native instrument plans")
        self.rank = rank

    def close(self) -> None:
        if self.context:
            self.library.sahay_destroy(self.context)
            self.context = None

    def marginal(self, bound: Bound, syndrome: Any, measured: int) -> float:
        return float(
            self.library.sahay_marginal(
                self.context, bound.coefficients, bound.values, syndrome, measured
            )
        )

    def amplitudes(self, bound: Bound, syndrome: Any) -> Any:
        result = np.empty(2, dtype=complex)
        self.library.sahay_amplitudes(
            self.context, bound.coefficients, bound.values, syndrome, result
        )
        return result
