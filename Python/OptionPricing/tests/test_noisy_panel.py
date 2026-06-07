import csv
import math

import numpy as np

from ImpliedVolatility.black_iv import implied_vol_black76
from ImpliedVolatility.black_price import black76_price, black76_vega
from OptionPricing.noisy_panel import (
    _apply_price_mechanics,
    default_noise_config,
    generate_noisy_panel_file,
    generate_noisy_panel_rows,
    read_table,
    write_noise_config,
    write_table,
)
from Scripts.validate_noisy_panels import validate_noisy_panels


def _clean_rows(sample_id=0):
    rows = []
    for week_index, spot in enumerate((100.0, 101.0)):
        for tau in (0.25, 0.5):
            forward = spot * math.exp(0.02 * tau)
            discount = math.exp(-0.02 * tau)
            for lm in (-0.05, 0.0, 0.05):
                strike = forward * math.exp(lm)
                option_type = "put" if lm < 0.0 else "call"
                vol = 0.20 + 0.02 * tau + 0.01 * abs(lm)
                price = black76_price(
                    forward=forward,
                    strike=strike,
                    tau=tau,
                    vol=vol,
                    discount_factor=discount,
                    option_type=option_type,
                )
                rows.append(
                    {
                        "run_id": "test_run",
                        "sample_id": sample_id,
                        "week_index": week_index,
                        "t": week_index / 52.0,
                        "S": spot,
                        "logS": math.log(spot),
                        "V": 0.04,
                        "r": 0.02,
                        "q": 0.0,
                        "maturity_years": tau,
                        "expiry_time": week_index / 52.0 + tau,
                        "forward": forward,
                        "log_moneyness": lm,
                        "strike": strike,
                        "option_type": option_type,
                        "is_otm": True,
                        "pricing_method": "COS",
                        "model_price": price,
                        "model_iv": implied_vol_black76(
                            price=price,
                            forward=forward,
                            strike=strike,
                            tau=tau,
                            discount_factor=discount,
                            option_type=option_type,
                        ),
                        "model_vega": black76_vega(
                            forward=forward,
                            strike=strike,
                            tau=tau,
                            vol=vol,
                            discount_factor=discount,
                        ),
                        "iv_method": "lets_be_rational",
                    }
                )
    return rows


def test_low_iid_noise_is_deterministic():
    config = default_noise_config()
    rows = _clean_rows()
    first, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=123, config=config)
    second, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=123, config=config)
    assert [row["raw_contaminated_iv"] for row in first] == [row["raw_contaminated_iv"] for row in second]
    assert all(row["noise_scenario"] == "low_iid" for row in first)
    assert all(row["observed_iv"] >= config.sigma_min for row in first)


def test_spatial_corr_noise_has_contract_level_draws():
    config = default_noise_config()
    rows, _ = generate_noisy_panel_rows(_clean_rows(), scenario="spatial_corr", seed=456, config=config)
    draws = np.array([row["noise_draw"] for row in rows])
    assert np.std(draws) > 0.0
    assert all(row["noise_scenario"] == "spatial_corr" for row in rows)


def test_persistent_factor_writes_factor_file_and_validation_passes(tmp_path):
    config = default_noise_config()
    run_root = tmp_path / "run"
    clean_path = write_table(_clean_rows(), run_root / "panels_clean" / "sample_000")
    results = []
    for scenario in config.enabled_scenarios():
        results.append(
            generate_noisy_panel_file(
                clean_panel_path=clean_path,
                run_root=run_root,
                sample_id=0,
                scenario=scenario,
                config=config,
            )
        )
    write_noise_config(run_root, config)
    assert {result.status for result in results} == {"ok"}
    factor_file = run_root / "noise_factors" / "persistent_factor" / "sample_000.csv"
    assert factor_file.exists()
    with factor_file.open(newline="") as fh:
        assert len(list(csv.DictReader(fh))) == 2
    validate_noisy_panels(run_root=run_root)


def test_tick_rounding_and_capping_are_recorded():
    config = default_noise_config()
    rows = _clean_rows()
    noisy, _ = generate_noisy_panel_rows(rows, scenario="low_iid", seed=1, config=config)
    assert all(abs(row["price_after_rounding"] / config.tick_size - round(row["price_after_rounding"] / config.tick_size)) < 1e-9 for row in noisy)

    _, _, was_capped, cap_direction = _apply_price_mechanics(
        raw_price=10_000.0,
        S=100.0,
        K=100.0,
        tau=0.25,
        r=0.02,
        q=0.0,
        option_type="call",
        tick_size=config.tick_size,
        price_epsilon=config.price_epsilon,
    )
    assert was_capped
    assert cap_direction == "upper"


def test_skip_existing_avoids_recomputation(tmp_path):
    config = default_noise_config()
    run_root = tmp_path / "run"
    clean_path = write_table(_clean_rows(), run_root / "panels_clean" / "sample_000")
    first = generate_noisy_panel_file(
        clean_panel_path=clean_path,
        run_root=run_root,
        sample_id=0,
        scenario="low_iid",
        config=config,
    )
    second = generate_noisy_panel_file(
        clean_panel_path=clean_path,
        run_root=run_root,
        sample_id=0,
        scenario="low_iid",
        config=config,
        skip_existing=True,
    )
    assert first.status == "ok"
    assert second.status == "skipped"
    assert len(read_table(first.output_observed_panel)) == len(_clean_rows())
