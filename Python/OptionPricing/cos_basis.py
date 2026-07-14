from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np

BASIS_CONVENTION = "fixed_effective_width"
PRICING_IMPLEMENTATION = "heston_cos_fixed_basis_v1"


def _canonical_json(payload: dict[str, Any]) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


@dataclass(frozen=True)
class FixedCosBasisConfig:
    """Shared maturity-specific COS widths and role-specific term counts."""

    maturities: tuple[float, ...]
    effective_widths: tuple[float, ...]
    generation_n_cos: int
    estimation_n_cos: int
    maturity_tolerance: float = 1e-10
    basis_convention: str = BASIS_CONVENTION
    pricing_implementation: str = PRICING_IMPLEMENTATION

    def validate(self) -> None:
        maturities = np.asarray(self.maturities, dtype=float)
        effective_widths = np.asarray(self.effective_widths, dtype=float)
        if self.basis_convention != BASIS_CONVENTION:
            raise ValueError(f"unsupported COS basis convention: {self.basis_convention!r}")
        if self.pricing_implementation != PRICING_IMPLEMENTATION:
            raise ValueError(f"unsupported COS pricing implementation: {self.pricing_implementation!r}")
        if maturities.ndim != 1 or effective_widths.ndim != 1 or maturities.size == 0:
            raise ValueError("maturities and effective_widths must be non-empty aligned 1D sequences")
        if maturities.shape != effective_widths.shape:
            raise ValueError("maturities and effective_widths must have a one-to-one match")
        if not np.all(np.isfinite(maturities)) or np.any(maturities <= 0.0):
            raise ValueError("fixed COS maturities must be finite and strictly positive")
        if not np.all(np.isfinite(effective_widths)) or np.any(effective_widths <= 0.0):
            raise ValueError("fixed COS effective widths must be finite and strictly positive")
        if self.generation_n_cos <= 1 or self.estimation_n_cos <= 1:
            raise ValueError("generation_n_cos and estimation_n_cos must both be > 1")
        if not np.isfinite(self.maturity_tolerance) or self.maturity_tolerance <= 0.0:
            raise ValueError("maturity_tolerance must be finite and strictly positive")
        ordered_maturities = np.sort(maturities)
        if ordered_maturities.size > 1 and np.any(np.diff(ordered_maturities) <= self.maturity_tolerance):
            raise ValueError("fixed COS maturity entries must be unique within maturity_tolerance")

    def width_for_maturity(self, maturity: float) -> float:
        self.validate()
        distances = np.abs(np.asarray(self.maturities, dtype=float) - float(maturity))
        matches = np.flatnonzero(distances <= self.maturity_tolerance)
        if matches.size == 0:
            raise ValueError(f"maturity {maturity!r} is not present in the fixed COS basis")
        if matches.size > 1:
            raise ValueError(f"maturity {maturity!r} matches multiple fixed COS basis entries")
        return float(self.effective_widths[int(matches[0])])

    def validate_requested_maturities(self, maturities: Iterable[float]) -> None:
        requested_maturities = np.asarray(tuple(maturities), dtype=float).reshape(-1)
        if requested_maturities.size == 0:
            raise ValueError("requested maturities must be non-empty")
        for maturity in requested_maturities:
            self.width_for_maturity(float(maturity))

    def to_dict(self) -> dict[str, Any]:
        self.validate()
        return {
            "basis_convention": self.basis_convention,
            "pricing_implementation": self.pricing_implementation,
            "maturities": [float(maturity) for maturity in self.maturities],
            "effective_widths": [float(width) for width in self.effective_widths],
            "generation_n_cos": int(self.generation_n_cos),
            "estimation_n_cos": int(self.estimation_n_cos),
            "maturity_tolerance": float(self.maturity_tolerance),
        }

    def stable_hash(self) -> str:
        return hashlib.sha256(_canonical_json(self.to_dict()).encode("utf-8")).hexdigest()

    def generation_metadata(self, *, resolved_basis_path: str | None = None) -> dict[str, Any]:
        metadata = {
            "basis_convention": self.basis_convention,
            "pricing_implementation": self.pricing_implementation,
            "maturities": [float(maturity) for maturity in self.maturities],
            "effective_widths": [float(width) for width in self.effective_widths],
            "generation_n_cos": int(self.generation_n_cos),
            "estimation_n_cos": int(self.estimation_n_cos),
            "maturity_tolerance": float(self.maturity_tolerance),
            "basis_config_hash": self.stable_hash(),
        }
        if resolved_basis_path is not None:
            metadata["resolved_basis_path"] = str(resolved_basis_path)
        return metadata

    def assert_panel_compatible(self, panel_metadata: dict[str, Any]) -> None:
        self.validate()
        if panel_metadata.get("basis_convention") != self.basis_convention:
            raise ValueError("panel and estimator COS basis conventions differ")
        if panel_metadata.get("pricing_implementation") != self.pricing_implementation:
            raise ValueError("panel and estimator COS pricing implementations differ")
        panel_maturities = tuple(float(value) for value in panel_metadata.get("maturities", ()))
        panel_widths = tuple(float(value) for value in panel_metadata.get("effective_widths", ()))
        if len(panel_maturities) != len(self.maturities) or len(panel_widths) != len(self.effective_widths):
            raise ValueError("panel and estimator COS maturity grids differ")
        if int(panel_metadata.get("generation_n_cos", 0)) != self.generation_n_cos:
            raise ValueError("panel generation_n_cos differs from the shared production basis")
        panel_tolerance = float(panel_metadata.get("maturity_tolerance", self.maturity_tolerance))
        comparison_tolerance = max(self.maturity_tolerance, panel_tolerance)
        for maturity, width in zip(panel_maturities, panel_widths):
            try:
                configured_width = self.width_for_maturity(maturity)
            except ValueError as error:
                raise ValueError("panel and estimator COS maturity grids differ") from error
            if not np.isclose(configured_width, width, rtol=0.0, atol=comparison_tolerance):
                raise ValueError(f"COS effective width mismatch at maturity {maturity}")
        self.validate_requested_maturities(panel_maturities)
        if panel_metadata.get("basis_config_hash") != self.stable_hash():
            raise ValueError("panel COS basis configuration hash differs from the estimator configuration")

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> FixedCosBasisConfig:
        config = cls(
            basis_convention=str(payload["basis_convention"]),
            pricing_implementation=str(payload["pricing_implementation"]),
            maturities=tuple(float(value) for value in payload["maturities"]),
            effective_widths=tuple(float(value) for value in payload["effective_widths"]),
            generation_n_cos=int(payload["generation_n_cos"]),
            estimation_n_cos=int(payload["estimation_n_cos"]),
            maturity_tolerance=float(payload.get("maturity_tolerance", 1e-10)),
        )
        config.validate()
        return config

    @classmethod
    def load(cls, file_path: str | Path) -> FixedCosBasisConfig:
        path = Path(file_path).expanduser().resolve()
        with path.open() as file_handle:
            payload = json.load(file_handle)
        if not isinstance(payload, dict):
            raise ValueError(f"COS basis file {path} must contain a JSON object")
        return cls.from_dict(payload)
