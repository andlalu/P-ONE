from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from OptionPricing.noisy_panel import (
    NOISE_SCENARIOS,
    NoiseSettings,
    clean_panel_file,
    generate_noisy_panel_file,
    observed_panel_file,
)
from Scripts.generation import generate_one_sample, load_config, set_thread_env, with_overrides


SCENARIO_LABELS = {
    "low_iid": "Design A: low i.i.d.",
    "spatial_corr": "Design B: spatial corr.",
    "persistent_factor": "Design C: persistent factor",
}

MATURITY_LABELS = {
    1.0 / 12.0: "1 month",
    0.25: "3 months",
    0.5: "6 months",
}

LINE_STYLES = {
    0: (0.00, None),
    1: (0.25, (5.0, 3.0)),
    2: (0.45, (1.0, 2.0)),
}

TIME_WINDOW_START = 100
TIME_WINDOW_LENGTH = 200


class PdfFigure:
    """Small single-page vector PDF helper for greyscale report figures."""

    def __init__(self, width_in: float, height_in: float) -> None:
        self.width = width_in * 72.0
        self.height = height_in * 72.0
        self.commands: list[str] = []

    @staticmethod
    def _num(value: float) -> str:
        return f"{value:.4f}".rstrip("0").rstrip(".")

    @staticmethod
    def _escape(text: str) -> str:
        return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    @staticmethod
    def _text_width(text: str, size: float, *, symbol: bool = False) -> float:
        width_factor = 0.55 if symbol else 0.50
        return width_factor * size * len(text)

    def rect(
        self,
        x: float,
        y: float,
        w: float,
        h: float,
        *,
        fill_gray: float | None = None,
        stroke_gray: float | None = 0.0,
        line_width: float = 0.4,
    ) -> None:
        self.commands.append("q")
        if fill_gray is not None:
            self.commands.append(f"{self._num(fill_gray)} g")
        if stroke_gray is not None:
            self.commands.append(f"{self._num(stroke_gray)} G {self._num(line_width)} w")
        self.commands.append(f"{self._num(x)} {self._num(y)} {self._num(w)} {self._num(h)} re")
        if fill_gray is not None and stroke_gray is not None:
            self.commands.append("B")
        elif fill_gray is not None:
            self.commands.append("f")
        elif stroke_gray is not None:
            self.commands.append("S")
        self.commands.append("Q")

    def line(
        self,
        x1: float,
        y1: float,
        x2: float,
        y2: float,
        *,
        gray: float = 0.0,
        width: float = 0.8,
        dash: tuple[float, ...] | None = None,
    ) -> None:
        self.commands.append("q")
        self.commands.append(f"{self._num(gray)} G {self._num(width)} w")
        self.commands.append("[] 0 d" if dash is None else f"[{' '.join(self._num(x) for x in dash)}] 0 d")
        self.commands.append(f"{self._num(x1)} {self._num(y1)} m {self._num(x2)} {self._num(y2)} l S")
        self.commands.append("Q")

    def polyline(
        self,
        points: list[tuple[float, float]],
        *,
        gray: float = 0.0,
        width: float = 0.9,
        dash: tuple[float, ...] | None = None,
    ) -> None:
        if len(points) < 2:
            return
        self.commands.append("q")
        self.commands.append(f"{self._num(gray)} G {self._num(width)} w")
        self.commands.append("[] 0 d" if dash is None else f"[{' '.join(self._num(x) for x in dash)}] 0 d")
        first = points[0]
        path = [f"{self._num(first[0])} {self._num(first[1])} m"]
        for x, y in points[1:]:
            path.append(f"{self._num(x)} {self._num(y)} l")
        self.commands.append(" ".join(path) + " S")
        self.commands.append("Q")

    def text(
        self,
        x: float,
        y: float,
        text: str,
        *,
        size: float = 8.0,
        gray: float = 0.0,
        align: str = "left",
        bold: bool = False,
        font_name: str | None = None,
    ) -> None:
        approx_width = self._text_width(text, size, symbol=font_name == "F3")
        if align == "center":
            x -= approx_width / 2.0
        elif align == "right":
            x -= approx_width
        font = font_name if font_name is not None else ("F2" if bold else "F1")
        self.commands.append("q")
        self.commands.append(f"{self._num(gray)} g")
        self.commands.append(f"BT /{font} {self._num(size)} Tf {self._num(x)} {self._num(y)} Td ({self._escape(text)}) Tj ET")
        self.commands.append("Q")

    @staticmethod
    def _math_chunks(text: str) -> list[tuple[str, str]]:
        macros = {"\\tau": "t", "\\Delta": "D"}
        chunks: list[tuple[str, str]] = []
        idx = 0
        while idx < len(text):
            matched = False
            for macro, symbol_char in macros.items():
                if text.startswith(macro, idx):
                    chunks.append(("symbol", symbol_char))
                    idx += len(macro)
                    matched = True
                    break
            if matched:
                continue
            next_idx = idx + 1
            while next_idx < len(text) and not any(text.startswith(macro, next_idx) for macro in macros):
                next_idx += 1
            chunks.append(("plain", text[idx:next_idx]))
            idx = next_idx
        return chunks

    def math_text(
        self,
        x: float,
        y: float,
        text: str,
        *,
        size: float = 8.0,
        gray: float = 0.0,
        align: str = "left",
        bold: bool = False,
    ) -> None:
        chunks = self._math_chunks(text)
        widths = [
            self._text_width(chunk, size, symbol=kind == "symbol")
            for kind, chunk in chunks
        ]
        total_width = sum(widths)
        if align == "center":
            x -= total_width / 2.0
        elif align == "right":
            x -= total_width
        for (kind, chunk), width in zip(chunks, widths):
            self.text(
                x,
                y,
                chunk,
                size=size,
                gray=gray,
                align="left",
                bold=bold,
                font_name="F3" if kind == "symbol" else None,
            )
            x += width

    def save(self, path: Path) -> None:
        content = ("\n".join(self.commands) + "\n").encode("latin-1", errors="replace")
        objects = [
            b"<< /Type /Catalog /Pages 2 0 R >>",
            b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
            (
                f"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {self._num(self.width)} {self._num(self.height)}] "
                "/Resources << /Font << /F1 4 0 R /F2 5 0 R /F3 6 0 R >> >> /Contents 7 0 R >>"
            ).encode("ascii"),
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold >>",
            b"<< /Type /Font /Subtype /Type1 /BaseFont /Symbol >>",
            f"<< /Length {len(content)} >>\nstream\n".encode("ascii") + content + b"endstream",
        ]
        out = bytearray(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n")
        offsets = [0]
        for idx, obj in enumerate(objects, start=1):
            offsets.append(len(out))
            out.extend(f"{idx} 0 obj\n".encode("ascii"))
            out.extend(obj)
            out.extend(b"\nendobj\n")
        xref_offset = len(out)
        out.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
        out.extend(b"0000000000 65535 f \n")
        for offset in offsets[1:]:
            out.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
        out.extend(
            (
                f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\n"
                f"startxref\n{xref_offset}\n%%EOF\n"
            ).encode("ascii")
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(bytes(out))


def _read_frame(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.suffix == ".csv":
        return pd.read_csv(path)
    raise ValueError(f"unsupported table format: {path}")


def _find_panel(base_without_suffix: Path) -> Path | None:
    for suffix in (".parquet", ".csv"):
        candidate = base_without_suffix.with_suffix(suffix)
        if candidate.exists():
            return candidate
    return None


def _required_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for name in candidates:
        if name in frame.columns:
            return name
    raise ValueError(f"missing {label} column; tried {', '.join(candidates)}")


def _panel_file_for_scenario(run_root: Path, scenario: str, sample_id: int) -> Path | None:
    return _find_panel(run_root / "panels_observed" / scenario / f"sample_{sample_id:03d}")


def _ensure_panels(
    *,
    config_path: Path,
    sample_id: int,
    generate_if_missing: bool,
) -> tuple[Path, dict[str, Path], Path | None, list[str], NoiseSettings]:
    config = load_config(config_path)
    if config.noise is None:
        raise ValueError("configuration does not contain noise scenarios")
    run_root = Path(config.output_root)
    generated: list[str] = []

    try:
        clean_path = clean_panel_file(run_root, sample_id)
    except FileNotFoundError:
        if not generate_if_missing:
            raise
        sample_config = with_overrides(
            config,
            n_samples=max(config.n_samples, sample_id + 1),
            workers=1,
            skip_existing=False,
        )
        result = generate_one_sample(sample_id, sample_config)
        if result.status == "error":
            raise RuntimeError(result.error)
        clean_path = Path(result.panel_file)
        generated.append("clean")
        generated.extend(config.noise.scenario_names())

    scenario_paths: dict[str, Path] = {}
    for scenario in NOISE_SCENARIOS:
        if scenario not in config.noise.scenario_names():
            continue
        existing = _panel_file_for_scenario(run_root, scenario, sample_id)
        factor_file = run_root / "noise_factors" / scenario / f"sample_{sample_id:03d}.csv"
        needs_factor = scenario == "persistent_factor" and not factor_file.exists()
        if existing is None or needs_factor:
            if not generate_if_missing:
                raise FileNotFoundError(f"missing observed panel for {scenario} sample {sample_id:03d}")
            result = generate_noisy_panel_file(
                clean_panel_path=clean_path,
                run_root=run_root,
                sample_id=sample_id,
                scenario=scenario,
                config=config.noise,
                panel_format=config.panel_format,
                skip_existing=False,
            )
            if result.status == "error":
                raise RuntimeError(result.error_message)
            existing = Path(result.output_observed_panel)
            generated.append(scenario)
        scenario_paths[scenario] = existing

    factor_path = run_root / "noise_factors" / "persistent_factor" / f"sample_{sample_id:03d}.csv"
    return clean_path, scenario_paths, factor_path if factor_path.exists() else None, generated, config.noise


def _date_noise_quality(observed: dict[str, pd.DataFrame], week_index: int) -> tuple[int, float]:
    capped = 0
    max_abs_delta_bp = 0.0
    for frame in observed.values():
        local = frame[frame["week_index"].astype(int) == int(week_index)]
        if local.empty:
            return 10**9, float("inf")
        delta = (local["observed_iv"].astype(float) - local["model_iv"].astype(float)) * 10000.0
        max_abs_delta_bp = max(max_abs_delta_bp, float(delta.abs().max()))
        if "was_price_capped" in local.columns:
            capped += int(local["was_price_capped"].astype(bool).sum())
    return capped, max_abs_delta_bp


def _select_representative_date(clean: pd.DataFrame, observed: dict[str, pd.DataFrame], requested: int) -> int:
    week_col = _required_column(clean, ("week_index",), "week index")
    weeks = sorted(int(x) for x in clean[week_col].drop_duplicates())
    if requested not in weeks:
        return weeks[len(weeks) // 2]

    capped, max_abs = _date_noise_quality(observed, requested)
    if capped == 0 and max_abs <= 150.0:
        return requested

    interior = weeks[5:-5] if len(weeks) > 10 else weeks
    candidates = []
    for week in interior:
        local_capped, local_max_abs = _date_noise_quality(observed, week)
        if local_capped == 0 and local_max_abs <= 150.0:
            candidates.append((abs(week - requested), week, local_max_abs))
    if candidates:
        return min(candidates)[1]
    return requested


def _date_slice(frame: pd.DataFrame, week_index: int) -> pd.DataFrame:
    out = frame[frame["week_index"].astype(int) == int(week_index)].copy()
    if out.empty:
        raise ValueError(f"date index {week_index} is unavailable in panel")
    return out.sort_values(["maturity_years", "log_moneyness"]).reset_index(drop=True)


def _grid_values(frame: pd.DataFrame, value_col: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    maturities = np.array(sorted(float(x) for x in frame["maturity_years"].drop_duplicates()), dtype=float)
    moneyness = np.array(sorted(float(x) for x in frame["log_moneyness"].drop_duplicates()), dtype=float)
    matrix = np.empty((len(maturities), len(moneyness)), dtype=float)
    for i, tau in enumerate(maturities):
        for j, lm in enumerate(moneyness):
            match = frame[
                np.isclose(frame["maturity_years"].astype(float), tau)
                & np.isclose(frame["log_moneyness"].astype(float), lm)
            ]
            if match.empty:
                raise ValueError(f"panel is missing maturity={tau} log_moneyness={lm}")
            matrix[i, j] = float(match.iloc[0][value_col])
    return maturities, moneyness, matrix


def _maturity_label(tau: float) -> str:
    for key, label in MATURITY_LABELS.items():
        if abs(tau - key) < 1.0e-8:
            return label
    return f"{tau:.3f}y"


def _nice_ticks(vmin: float, vmax: float, count: int = 4) -> list[float]:
    if not math.isfinite(vmin) or not math.isfinite(vmax):
        return [0.0]
    if abs(vmax - vmin) < 1.0e-12:
        margin = max(abs(vmin), 1.0) * 0.05
        vmin -= margin
        vmax += margin
    raw = (vmax - vmin) / max(count - 1, 1)
    power = 10.0 ** math.floor(math.log10(abs(raw)))
    step = min((1.0, 2.0, 2.5, 5.0, 10.0), key=lambda item: abs(item * power - raw)) * power
    start = math.floor(vmin / step) * step
    end = math.ceil(vmax / step) * step
    ticks = []
    value = start
    while value <= end + 0.5 * step:
        ticks.append(0.0 if abs(value) < 1.0e-12 else value)
        value += step
    return ticks


def _fmt_decimal(value: float) -> str:
    return f"{value:.3f}".rstrip("0").rstrip(".")


def _fmt_bp(value: float) -> str:
    return f"{value:.1f}"


def _draw_axes(
    pdf: PdfFigure,
    rect: tuple[float, float, float, float],
    *,
    xmin: float,
    xmax: float,
    ymin: float,
    ymax: float,
    title: str,
    xlabel: str = "",
    ylabel: str = "",
    x_ticks: list[float] | None = None,
    y_ticks: list[float] | None = None,
) -> tuple[Callable[[float], float], Callable[[float], float]]:
    x0, y0, w, h = rect
    pdf.rect(x0, y0, w, h, fill_gray=None, stroke_gray=0.0, line_width=0.5)
    x_ticks = x_ticks if x_ticks is not None else _nice_ticks(xmin, xmax, 5)
    y_ticks = y_ticks if y_ticks is not None else _nice_ticks(ymin, ymax, 4)

    def sx(value: float) -> float:
        return x0 + (value - xmin) / (xmax - xmin) * w

    def sy(value: float) -> float:
        return y0 + (value - ymin) / (ymax - ymin) * h

    for value in x_ticks:
        if xmin - 1.0e-12 <= value <= xmax + 1.0e-12:
            x = sx(value)
            pdf.line(x, y0, x, y0 - 3.0, gray=0.0, width=0.4)
            pdf.text(x, y0 - 12.0, _fmt_decimal(value), size=6.5, align="center")
    for value in y_ticks:
        if ymin - 1.0e-12 <= value <= ymax + 1.0e-12:
            y = sy(value)
            pdf.line(x0 - 3.0, y, x0, y, gray=0.0, width=0.4)
            pdf.text(x0 - 5.0, y - 2.2, _fmt_decimal(value), size=6.5, align="right")
            pdf.line(x0, y, x0 + w, y, gray=0.88, width=0.25)
    pdf.text(x0 + w / 2.0, y0 + h + 9.0, title, size=8.5, align="center", bold=True)
    if xlabel:
        pdf.math_text(x0 + w / 2.0, y0 - 25.0, xlabel, size=7.2, align="center")
    if ylabel:
        pdf.math_text(x0 - 34.0, y0 + h / 2.0, ylabel, size=7.2, align="center")
    return sx, sy


def _draw_legend(
    pdf: PdfFigure,
    x: float,
    y: float,
    labels: list[str],
    *,
    width: float = 42.0,
) -> None:
    pdf.rect(x - 4.0, y - len(labels) * 10.0 + 3.0, width, len(labels) * 10.0 + 8.0, fill_gray=1.0, stroke_gray=0.82, line_width=0.25)
    for idx, label in enumerate(labels):
        gray, dash = LINE_STYLES[idx % len(LINE_STYLES)]
        yy = y - idx * 10.0
        pdf.line(x, yy + 3.0, x + 20.0, yy + 3.0, gray=gray, width=1.1, dash=dash)
        pdf.math_text(x + 24.0, yy, label, size=7.0)


def _draw_legend_horizontal(pdf: PdfFigure, x: float, y: float, labels: list[str]) -> None:
    offset = 0.0
    for idx, label in enumerate(labels):
        gray, dash = LINE_STYLES[idx % len(LINE_STYLES)]
        pdf.line(x + offset, y + 3.0, x + offset + 20.0, y + 3.0, gray=gray, width=1.1, dash=dash)
        pdf.text(x + offset + 24.0, y, label, size=7.0)
        offset += 88.0


def _line_panel(
    pdf: PdfFigure,
    rect: tuple[float, float, float, float],
    *,
    frame: pd.DataFrame,
    value_col: str,
    title: str,
    y_min: float,
    y_max: float,
    clean_reference: pd.DataFrame | None = None,
    show_xlabel: bool = False,
    show_ylabel: bool = False,
) -> None:
    maturities, moneyness, _ = _grid_values(frame, value_col)
    sx, sy = _draw_axes(
        pdf,
        rect,
        xmin=float(moneyness.min()),
        xmax=float(moneyness.max()),
        ymin=y_min,
        ymax=y_max,
        title=title,
        xlabel="log-forward moneyness" if show_xlabel else "",
        ylabel="implied vol." if show_ylabel else "",
        x_ticks=[float(x) for x in moneyness],
        y_ticks=_nice_ticks(y_min, y_max, 4),
    )
    if clean_reference is not None:
        for tau in maturities:
            ref = clean_reference[np.isclose(clean_reference["maturity_years"].astype(float), tau)].sort_values("log_moneyness")
            points = [(sx(float(row["log_moneyness"])), sy(float(row["model_iv"]))) for _, row in ref.iterrows()]
            pdf.polyline(points, gray=0.72, width=0.6, dash=(2.0, 2.0))

    for idx, tau in enumerate(maturities):
        part = frame[np.isclose(frame["maturity_years"].astype(float), tau)].sort_values("log_moneyness")
        gray, dash = LINE_STYLES[idx % len(LINE_STYLES)]
        points = [(sx(float(row["log_moneyness"])), sy(float(row[value_col]))) for _, row in part.iterrows()]
        pdf.polyline(points, gray=gray, width=1.0, dash=dash)


def _draw_heatmap(
    pdf: PdfFigure,
    rect: tuple[float, float, float, float],
    matrix: np.ndarray,
    *,
    x_labels: list[str],
    y_labels: list[str],
    title: str,
    gray_for_value: Callable[[float], float],
    value_label: Callable[[float], str] | None = None,
    xlabel: str = "",
    ylabel: str = "",
) -> None:
    x0, y0, w, h = rect
    rows, cols = matrix.shape
    cell_w = w / cols
    cell_h = h / rows
    pdf.text(x0 + w / 2.0, y0 + h + 10.0, title, size=8.5, align="center", bold=True)
    for i in range(rows):
        for j in range(cols):
            value = float(matrix[i, j])
            gray = max(0.0, min(1.0, gray_for_value(value)))
            yy = y0 + (rows - 1 - i) * cell_h
            xx = x0 + j * cell_w
            pdf.rect(xx, yy, cell_w, cell_h, fill_gray=gray, stroke_gray=0.18, line_width=0.25)
            if value_label is not None:
                text_gray = 1.0 if gray < 0.35 else 0.0
                pdf.text(xx + cell_w / 2.0, yy + cell_h / 2.0 - 2.2, value_label(value), size=6.2, gray=text_gray, align="center")
    for j, label in enumerate(x_labels):
        pdf.text(x0 + (j + 0.5) * cell_w, y0 - 11.0, label, size=6.5, align="center")
    for i, label in enumerate(y_labels):
        yy = y0 + (rows - 1 - i + 0.5) * cell_h
        pdf.text(x0 - 5.0, yy - 2.2, label, size=6.5, align="right")
    if xlabel:
        pdf.math_text(x0 + w / 2.0, y0 - 24.0, xlabel, size=7.2, align="center")
    if ylabel:
        pdf.math_text(x0 - 34.0, y0 + h / 2.0, ylabel, size=7.2, align="center")


def _draw_colorbar(
    pdf: PdfFigure,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    min_label: str,
    mid_label: str | None,
    max_label: str,
    title: str,
    gray_fn: Callable[[float], float],
) -> None:
    steps = 40
    for idx in range(steps):
        t0 = idx / steps
        gray = gray_fn(t0)
        pdf.rect(x + t0 * w, y, w / steps + 0.2, h, fill_gray=gray, stroke_gray=None)
    pdf.rect(x, y, w, h, fill_gray=None, stroke_gray=0.0, line_width=0.25)
    pdf.text(x, y - 10.0, min_label, size=6.5, align="left")
    if mid_label is not None:
        pdf.text(x + w / 2.0, y - 10.0, mid_label, size=6.5, align="center")
    pdf.text(x + w, y - 10.0, max_label, size=6.5, align="right")
    pdf.math_text(x + w / 2.0, y + h + 4.0, title, size=6.8, align="center")


def _scale_gray(min_value: float, max_value: float) -> Callable[[float], float]:
    denom = max(max_value - min_value, 1.0e-12)

    def mapper(value: float) -> float:
        z = (value - min_value) / denom
        return 0.92 - 0.76 * max(0.0, min(1.0, z))

    return mapper


def _diverging_gray(vlim: float) -> Callable[[float], float]:
    limit = max(abs(vlim), 1.0e-12)

    def mapper(value: float) -> float:
        z = (value + limit) / (2.0 * limit)
        return 0.90 - 0.78 * max(0.0, min(1.0, z))

    return mapper


def _make_surface_figure(
    output: Path,
    *,
    clean_date: pd.DataFrame,
    observed_dates: dict[str, pd.DataFrame],
) -> None:
    values = [clean_date["model_iv"].astype(float).to_numpy()]
    for frame in observed_dates.values():
        values.append(frame["observed_iv"].astype(float).to_numpy())
    all_values = np.concatenate(values)
    y_min = float(np.nanmin(all_values))
    y_max = float(np.nanmax(all_values))
    pad = max(0.0025, 0.08 * (y_max - y_min))
    y_min -= pad
    y_max += pad

    pdf = PdfFigure(7.2, 5.2)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Clean and contaminated implied-volatility panels", size=10.0, align="center", bold=True)
    left = 72.0
    right = 18.0
    bottom = 66.0
    top = 62.0
    gap_x = 34.0
    gap_y = 42.0
    panel_w = (pdf.width - left - right - gap_x) / 2.0
    panel_h = (pdf.height - bottom - top - gap_y) / 2.0
    rects = [
        (left, bottom + panel_h + gap_y, panel_w, panel_h),
        (left + panel_w + gap_x, bottom + panel_h + gap_y, panel_w, panel_h),
        (left, bottom, panel_w, panel_h),
        (left + panel_w + gap_x, bottom, panel_w, panel_h),
    ]
    _line_panel(
        pdf,
        rects[0],
        frame=clean_date,
        value_col="model_iv",
        title="Clean panel",
        y_min=y_min,
        y_max=y_max,
        clean_reference=None,
        show_ylabel=True,
    )
    for rect, scenario in zip(rects[1:], ("low_iid", "spatial_corr", "persistent_factor")):
        _line_panel(
            pdf,
            rect,
            frame=observed_dates[scenario],
            value_col="observed_iv",
            title=SCENARIO_LABELS[scenario],
            y_min=y_min,
            y_max=y_max,
            clean_reference=clean_date,
            show_xlabel=rect[1] == bottom,
            show_ylabel=rect[0] == left,
        )
    maturities = [float(x) for x in sorted(clean_date["maturity_years"].drop_duplicates())]
    _draw_legend_horizontal(pdf, 143.0, pdf.height - 38.0, [_maturity_label(tau) for tau in maturities])
    pdf.text(pdf.width - 18.0, 14.0, "Light grey curves in noisy panels: clean reference.", size=6.8, align="right")
    pdf.save(output)


def _make_difference_figure(
    output: Path,
    *,
    observed_dates: dict[str, pd.DataFrame],
) -> None:
    diff_mats: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    all_diffs: list[np.ndarray] = []
    for scenario, frame in observed_dates.items():
        local = frame.copy()
        local["delta_iv_bp"] = (local["observed_iv"].astype(float) - local["model_iv"].astype(float)) * 10000.0
        mats = _grid_values(local, "delta_iv_bp")
        diff_mats[scenario] = mats
        all_diffs.append(mats[2])
    vlim = float(np.nanmax(np.abs(np.concatenate([x.ravel() for x in all_diffs]))))
    mapper = _diverging_gray(vlim)

    pdf = PdfFigure(7.2, 4.8)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Signed implied-volatility perturbations", size=10.0, align="center", bold=True)
    left = 46.0
    right = 28.0
    bottom = 96.0
    gap = 24.0
    panel_w = (pdf.width - left - right - 2.0 * gap) / 3.0
    panel_h = 150.0
    for idx, scenario in enumerate(("low_iid", "spatial_corr", "persistent_factor")):
        maturities, moneyness, matrix = diff_mats[scenario]
        rect = (left + idx * (panel_w + gap), bottom, panel_w, panel_h)
        _draw_heatmap(
            pdf,
            rect,
            matrix,
            x_labels=[_fmt_decimal(x) for x in moneyness],
            y_labels=[_maturity_label(x).replace(" months", "m").replace(" month", "m") for x in maturities],
            title=SCENARIO_LABELS[scenario],
            gray_for_value=mapper,
            value_label=_fmt_bp,
            xlabel="m" if idx == 1 else "",
            ylabel="\\tau" if idx == 0 else "",
        )
    _draw_colorbar(
        pdf,
        170.0,
        44.0,
        180.0,
        9.0,
        min_label=f"-{vlim:.1f}",
        mid_label="0",
        max_label=f"{vlim:.1f}",
        title="\\Delta IV, basis points",
        gray_fn=lambda t: 0.90 - 0.78 * t,
    )
    pdf.save(output)


def _marginal_scale(log_moneyness: np.ndarray, tau: np.ndarray, config: Any) -> np.ndarray:
    return config.alpha_0 * (
        1.0 + config.alpha_m * np.abs(log_moneyness) + config.alpha_tau / np.sqrt(np.maximum(tau, config.tau_min))
    )


def _make_structure_figure(
    output: Path,
    *,
    clean: pd.DataFrame,
    noise_config: NoiseSettings,
) -> None:
    base = _date_slice(clean, int(clean["week_index"].iloc[0]))
    maturities = np.array(sorted(float(x) for x in base["maturity_years"].drop_duplicates()), dtype=float)
    moneyness = np.array(sorted(float(x) for x in base["log_moneyness"].drop_duplicates()), dtype=float)
    scale = np.empty((len(maturities), len(moneyness)), dtype=float)
    for i, tau in enumerate(maturities):
        for j, lm in enumerate(moneyness):
            scale[i, j] = float(
                _marginal_scale(
                    np.array([lm]),
                    np.array([tau]),
                    noise_config.scenarios["low_iid"],
                )[0]
            ) * 10000.0

    ordered = base.sort_values(["maturity_years", "log_moneyness"]).reset_index(drop=True)
    lm = ordered["log_moneyness"].astype(float).to_numpy()
    tau = ordered["maturity_years"].astype(float).to_numpy()
    spatial = noise_config.scenarios["spatial_corr"]
    distance = np.abs(lm[:, None] - lm[None, :]) / float(spatial["ell_m"])
    distance += np.abs(tau[:, None] - tau[None, :]) / float(spatial["ell_tau"])
    corr = np.exp(-distance)

    pdf = PdfFigure(7.2, 4.8)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, "Noise scale and spatial dependence structure", size=10.0, align="center", bold=True)
    scale_mapper = _scale_gray(float(scale.min()), float(scale.max()))
    corr_mapper = _scale_gray(0.0, 1.0)
    _draw_heatmap(
        pdf,
        (58.0, 118.0, 185.0, 135.0),
        scale,
        x_labels=[_fmt_decimal(x) for x in moneyness],
        y_labels=[_maturity_label(x).replace(" months", "m").replace(" month", "m") for x in maturities],
        title="Marginal noise scale",
        gray_for_value=scale_mapper,
        value_label=lambda x: f"{x:.1f}",
        xlabel="log-forward moneyness",
        ylabel="\\tau",
    )
    _draw_colorbar(
        pdf,
        88.0,
        72.0,
        120.0,
        9.0,
        min_label=f"{scale.min():.1f}",
        mid_label=None,
        max_label=f"{scale.max():.1f}",
        title="IV bp",
        gray_fn=lambda t: 0.92 - 0.76 * t,
    )
    _draw_heatmap(
        pdf,
        (306.0, 72.0, 180.0, 180.0),
        corr,
        x_labels=[str(i + 1) if i in (0, 4, 9, 14) else "" for i in range(corr.shape[1])],
        y_labels=[str(i + 1) if i in (0, 4, 9, 14) else "" for i in range(corr.shape[0])],
        title="Design B correlation matrix",
        gray_for_value=corr_mapper,
        value_label=None,
        xlabel="",
        ylabel="contract index",
    )
    _draw_colorbar(
        pdf,
        338.0,
        28.0,
        120.0,
        9.0,
        min_label="0",
        mid_label="0.5",
        max_label="1",
        title="correlation",
        gray_fn=lambda t: 0.92 - 0.76 * t,
    )
    pdf.save(output)


def _resolve_time_window(frame: pd.DataFrame, start: int, length: int) -> tuple[int, int]:
    weeks = sorted(int(x) for x in frame["week_index"].drop_duplicates())
    if not weeks:
        raise ValueError("cannot draw time-series figure from an empty panel")
    max_start = max(weeks[-1] - length + 1, weeks[0])
    window_start = min(max(start, weeks[0]), max_start)
    window_end = min(window_start + length - 1, weeks[-1])
    return int(window_start), int(window_end)


def _week_ticks(start: int, end: int) -> list[float]:
    step = 50
    first = int(math.ceil(start / step) * step)
    ticks = [float(start)]
    ticks.extend(float(value) for value in range(first, end + 1, step) if value != start)
    if end not in {int(value) for value in ticks}:
        ticks.append(float(end))
    return ticks


def _selected_contract_series(
    observed: pd.DataFrame,
    *,
    start: int,
    end: int,
) -> list[tuple[str, np.ndarray, np.ndarray]]:
    local = observed.copy()
    local["delta_iv_bp"] = (local["observed_iv"].astype(float) - local["model_iv"].astype(float)) * 10000.0
    local = local[(local["week_index"].astype(int) >= start) & (local["week_index"].astype(int) <= end)]
    contract_specs = [
        ("ATM short", 0.0, observed["maturity_years"].astype(float).min()),
        ("Low m, 3m", observed["log_moneyness"].astype(float).min(), 0.25),
        ("High m, 6m", observed["log_moneyness"].astype(float).max(), observed["maturity_years"].astype(float).max()),
    ]
    series: list[tuple[str, np.ndarray, np.ndarray]] = []
    for label, lm, tau in contract_specs:
        part = local[
            np.isclose(local["log_moneyness"].astype(float), lm)
            & np.isclose(local["maturity_years"].astype(float), tau)
        ].sort_values("week_index")
        if part.empty:
            raise ValueError(f"missing representative contract series for {label}")
        series.append((label, part["week_index"].astype(float).to_numpy(), part["delta_iv_bp"].astype(float).to_numpy()))
    return series


def _make_low_iid_timeseries_figure(
    output: Path,
    *,
    observed: pd.DataFrame,
    window_start: int = TIME_WINDOW_START,
    window_length: int = TIME_WINDOW_LENGTH,
) -> None:
    start, end = _resolve_time_window(observed, window_start, window_length)
    series = _selected_contract_series(observed, start=start, end=end)
    all_y = np.concatenate([item[2] for item in series])
    vlim = max(1.0, float(np.nanmax(np.abs(all_y))) * 1.08)

    pdf = PdfFigure(7.2, 3.55)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, f"Design A i.i.d. perturbations, weeks {start}-{end}", size=10.0, align="center", bold=True)
    rect = (58.0, 70.0, 417.0, 145.0)
    sx, sy = _draw_axes(
        pdf,
        rect,
        xmin=float(start),
        xmax=float(end),
        ymin=-vlim,
        ymax=vlim,
        title="Selected contract perturbations",
        xlabel="week index",
        ylabel="\\Delta IV bp",
        x_ticks=_week_ticks(start, end),
        y_ticks=_nice_ticks(-vlim, vlim, 5),
    )
    pdf.line(rect[0], sy(0.0), rect[0] + rect[2], sy(0.0), gray=0.55, width=0.35, dash=(2.0, 2.0))
    for idx, (_, x_values, y_values) in enumerate(series):
        gray, dash = LINE_STYLES[idx]
        points = [(sx(float(xx)), sy(float(yy))) for xx, yy in zip(x_values, y_values)]
        pdf.polyline(points, gray=gray, width=0.9, dash=dash)
    _draw_legend(pdf, 390.0, 192.0, [item[0] for item in series], width=88.0)
    pdf.save(output)


def _make_persistent_factor_figure(
    output: Path,
    *,
    observed: pd.DataFrame,
    factor_file: Path | None,
    window_start: int = TIME_WINDOW_START,
    window_length: int = TIME_WINDOW_LENGTH,
) -> None:
    start, end = _resolve_time_window(observed, window_start, window_length)
    pdf = PdfFigure(7.2, 5.0)
    pdf.text(pdf.width / 2.0, pdf.height - 18.0, f"Persistent-factor observation-error dynamics, weeks {start}-{end}", size=10.0, align="center", bold=True)
    top_rect = (55.0, 214.0, 420.0, 95.0)
    bottom_rect = (55.0, 58.0, 420.0, 105.0)

    if factor_file is not None:
        factors = pd.read_csv(factor_file)
        factors = factors[(factors["week_index"].astype(int) >= start) & (factors["week_index"].astype(int) <= end)]
        x = factors["week_index"].astype(float).to_numpy()
        series = [
            ("f0 level", factors["factor_0"].astype(float).to_numpy() * 10000.0),
            ("f1 slope m", factors["factor_1"].astype(float).to_numpy() * 10000.0),
            ("f2 slope \\tau", factors["factor_2"].astype(float).to_numpy() * 10000.0),
        ]
        all_y = np.concatenate([item[1] for item in series])
        pad = max(1.0, 0.08 * (float(all_y.max()) - float(all_y.min())))
        sx, sy = _draw_axes(
            pdf,
            top_rect,
            xmin=float(start),
            xmax=float(end),
            ymin=float(all_y.min() - pad),
            ymax=float(all_y.max() + pad),
            title="Realised common factors",
            xlabel="week index",
            ylabel="IV bp",
            x_ticks=_week_ticks(start, end),
            y_ticks=_nice_ticks(float(all_y.min() - pad), float(all_y.max() + pad), 4),
        )
        for idx, (_, y_values) in enumerate(series):
            gray, dash = LINE_STYLES[idx]
            points = [(sx(float(xx)), sy(float(yy))) for xx, yy in zip(x, y_values)]
            pdf.polyline(points, gray=gray, width=0.9, dash=dash)
        _draw_legend(pdf, 395.0, 294.0, [item[0] for item in series], width=88.0)
    else:
        pdf.text(top_rect[0] + top_rect[2] / 2.0, top_rect[1] + top_rect[3] / 2.0, "Persistent factor file unavailable", size=8.0, align="center")

    series2 = _selected_contract_series(observed, start=start, end=end)
    all_y2 = np.concatenate([item[2] for item in series2])
    pad2 = max(1.0, 0.08 * (float(all_y2.max()) - float(all_y2.min())))
    sx2, sy2 = _draw_axes(
        pdf,
        bottom_rect,
        xmin=float(start),
        xmax=float(end),
        ymin=float(all_y2.min() - pad2),
        ymax=float(all_y2.max() + pad2),
        title="Selected contract perturbations",
        xlabel="week index",
        ylabel="\\Delta IV bp",
        x_ticks=_week_ticks(start, end),
        y_ticks=_nice_ticks(float(all_y2.min() - pad2), float(all_y2.max() + pad2), 4),
    )
    pdf.line(bottom_rect[0], sy2(0.0), bottom_rect[0] + bottom_rect[2], sy2(0.0), gray=0.55, width=0.35, dash=(2.0, 2.0))
    for idx, (_, x_values, y_values) in enumerate(series2):
        gray, dash = LINE_STYLES[idx]
        points = [(sx2(float(xx)), sy2(float(yy))) for xx, yy in zip(x_values, y_values)]
        pdf.polyline(points, gray=gray, width=0.9, dash=dash)
    _draw_legend(pdf, 395.0, 140.0, [item[0] for item in series2], width=88.0)
    pdf.save(output)


def _write_manifest(
    path: Path,
    *,
    sample_id: int,
    requested_date_index: int,
    selected_date_index: int,
    clean_path: Path,
    scenario_paths: dict[str, Path],
    factor_file: Path | None,
    figures: list[Path],
    generated: list[str],
    time_window: tuple[int, int],
) -> None:
    payload = {
        "sample_id": sample_id,
        "requested_date_index": requested_date_index,
        "date_index": selected_date_index,
        "time_window": {
            "start_week_index": time_window[0],
            "end_week_index": time_window[1],
            "length": time_window[1] - time_window[0] + 1,
        },
        "scenario_files": {
            "clean": str(clean_path),
            **{scenario: str(scenario_paths[scenario]) for scenario in sorted(scenario_paths)},
        },
        "factor_file": str(factor_file) if factor_file is not None else None,
        "figures": [str(item) for item in figures],
        "generated_or_regenerated": generated,
        "created_by": "Python/Scripts/make_noise_visualisations.py",
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Create greyscale noise-panel visualisations for the June 2026 report.")
    parser.add_argument("--config", default="Python/Scripts/configs/heston_experiment_run_001.json")
    parser.add_argument("--sample-id", type=int, default=0)
    parser.add_argument("--date-index", type=int, default=100)
    parser.add_argument("--output-dir", default="Reports/18062026/figures/noise")
    parser.add_argument("--generate-if-missing", action="store_true")
    args = parser.parse_args()

    set_thread_env()
    clean_path, scenario_paths, factor_file, generated, noise_config = _ensure_panels(
        config_path=Path(args.config),
        sample_id=args.sample_id,
        generate_if_missing=args.generate_if_missing,
    )
    clean = _read_frame(clean_path)
    observed = {scenario: _read_frame(path) for scenario, path in scenario_paths.items()}
    selected_date = _select_representative_date(clean, observed, args.date_index)
    if selected_date != args.date_index:
        available_weeks = set(int(x) for x in clean["week_index"].drop_duplicates())
        if args.date_index not in available_weeks:
            print(f"requested date-index {args.date_index} unavailable; selected interior date-index {selected_date}")
        else:
            capped, max_abs = _date_noise_quality(observed, args.date_index)
            print(
                f"requested date-index {args.date_index} is cap-dominated "
                f"(capped={capped}, max_abs_delta_bp={max_abs:.1f}); selected representative date-index {selected_date}"
            )
    else:
        print(f"selected date-index {selected_date}")

    clean_date = _date_slice(clean, selected_date)
    observed_dates = {scenario: _date_slice(frame, selected_date) for scenario, frame in observed.items()}
    missing = sorted(set(NOISE_SCENARIOS) - set(observed_dates))
    if missing:
        raise ValueError(f"missing enabled noisy panels: {', '.join(missing)}")

    output_dir = Path(args.output_dir)
    sample_tag = f"sample{args.sample_id:03d}"
    date_tag = f"date{selected_date:03d}"
    time_window = _resolve_time_window(observed["low_iid"], TIME_WINDOW_START, TIME_WINDOW_LENGTH)
    figures = [
        output_dir / f"noise_surface_comparison_{sample_tag}_{date_tag}.pdf",
        output_dir / f"noise_difference_comparison_{sample_tag}_{date_tag}.pdf",
        output_dir / "noise_design_cross_section_structure.pdf",
        output_dir / f"noise_low_iid_timeseries_{sample_tag}.pdf",
        output_dir / f"noise_persistent_factor_timeseries_{sample_tag}.pdf",
    ]
    _make_surface_figure(figures[0], clean_date=clean_date, observed_dates=observed_dates)
    _make_difference_figure(figures[1], observed_dates=observed_dates)
    _make_structure_figure(figures[2], clean=clean, noise_config=noise_config)
    _make_low_iid_timeseries_figure(figures[3], observed=observed["low_iid"])
    _make_persistent_factor_figure(figures[4], observed=observed["persistent_factor"], factor_file=factor_file)
    _write_manifest(
        output_dir / "noise_visualisation_manifest.json",
        sample_id=args.sample_id,
        requested_date_index=args.date_index,
        selected_date_index=selected_date,
        clean_path=clean_path,
        scenario_paths=scenario_paths,
        factor_file=factor_file,
        figures=figures,
        generated=generated,
        time_window=time_window,
    )

    print(f"sample-id {args.sample_id}")
    print(f"time window: week-index {time_window[0]} to {time_window[1]}")
    print(f"clean panel: {clean_path}")
    for scenario in NOISE_SCENARIOS:
        print(f"{scenario} panel: {scenario_paths[scenario]}")
    print(f"factor file: {factor_file if factor_file is not None else 'not available'}")
    for figure in figures:
        print(f"wrote {figure}")
    print(f"manifest: {output_dir / 'noise_visualisation_manifest.json'}")
    if generated:
        print(f"generated or regenerated: {', '.join(generated)}")
    else:
        print("loaded existing panels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
