#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class MetricPolicy:
    name: str
    description: str
    max_regression_pct: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--policy", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--current", required=True)
    parser.add_argument("--baseline")
    parser.add_argument("--report-json")
    parser.add_argument("--strict-metrics", action="store_true")
    return parser.parse_args()


def load_policy(policy_path: Path, profile_name: str) -> tuple[Path | None, float, dict[str, MetricPolicy]]:
    payload = json.loads(policy_path.read_text(encoding="utf-8"))
    profiles = payload.get("profiles", {})
    profile = profiles.get(profile_name)
    if profile is None:
        raise SystemExit(f"profilo benchmark sconosciuto: {profile_name}")

    default_threshold = float(payload.get("default_threshold_pct", 15.0))
    baseline_value = profile.get("baseline")
    baseline_path = (policy_path.parent.parent.parent / baseline_value).resolve() if baseline_value else None

    metrics: dict[str, MetricPolicy] = {}
    for name, data in profile.get("metrics", {}).items():
        metrics[name] = MetricPolicy(
            name=name,
            description=str(data.get("description", "")),
            max_regression_pct=float(data.get("max_regression_pct", default_threshold)),
        )
    return baseline_path, default_threshold, metrics


def read_metrics(path: Path) -> dict[str, float]:
    metrics: dict[str, float] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        parts = raw_line.strip().split()
        if len(parts) != 2:
            continue
        name, value = parts
        metrics[name] = float(value)
    return metrics


def compare_metrics(
    baseline: dict[str, float],
    current: dict[str, float],
    metric_policies: dict[str, MetricPolicy],
    default_threshold_pct: float,
    strict_metrics: bool,
) -> dict[str, Any]:
    missing_metrics = sorted(name for name in baseline if name not in current)
    unexpected_metrics = sorted(name for name in current if name not in baseline)
    results: list[dict[str, Any]] = []
    has_regression = False

    metric_names = sorted(set(baseline.keys()) | set(current.keys()))
    for name in metric_names:
        baseline_value = baseline.get(name)
        current_value = current.get(name)
        metric_policy = metric_policies.get(
            name,
            MetricPolicy(name=name, description="", max_regression_pct=default_threshold_pct),
        )

        status = "ok"
        allowed_value = None
        delta_pct = None
        if baseline_value is None:
            status = "new"
        elif current_value is None:
            status = "missing"
            has_regression = True
        else:
            allowed_value = baseline_value * (1.0 + metric_policy.max_regression_pct / 100.0)
            if baseline_value == 0.0:
                delta_pct = math.inf if current_value > 0.0 else 0.0
            else:
                delta_pct = ((current_value / baseline_value) - 1.0) * 100.0
            if current_value > allowed_value:
                status = "regression"
                has_regression = True

        results.append(
            {
                "name": name,
                "description": metric_policy.description,
                "baseline_ms": baseline_value,
                "current_ms": current_value,
                "allowed_ms": allowed_value,
                "delta_pct": delta_pct,
                "threshold_pct": metric_policy.max_regression_pct,
                "status": status,
            }
        )

    if strict_metrics and unexpected_metrics:
        has_regression = True

    return {
        "status": "failed" if has_regression else "ok",
        "strict_metrics": strict_metrics,
        "missing_metrics": missing_metrics,
        "unexpected_metrics": unexpected_metrics,
        "metrics": results,
    }


def emit_text_report(report: dict[str, Any]) -> None:
    for metric in report["metrics"]:
        baseline_value = metric["baseline_ms"]
        current_value = metric["current_ms"]
        allowed_value = metric["allowed_ms"]
        delta_pct = metric["delta_pct"]

        baseline_text = "-" if baseline_value is None else f"{baseline_value:.3f}"
        current_text = "-" if current_value is None else f"{current_value:.3f}"
        allowed_text = "-" if allowed_value is None else f"{allowed_value:.3f}"
        if delta_pct is None:
            delta_text = "-"
        elif math.isinf(delta_pct):
            delta_text = "inf"
        else:
            delta_text = f"{delta_pct:.2f}"

        print(
            "benchmark_report "
            f"name={metric['name']} "
            f"status={metric['status']} "
            f"baseline_ms={baseline_text} "
            f"current_ms={current_text} "
            f"allowed_ms={allowed_text} "
            f"delta_pct={delta_text} "
            f"threshold_pct={metric['threshold_pct']:.2f}"
        )

    if report["missing_metrics"]:
        print("benchmark_report missing_metrics=" + ",".join(report["missing_metrics"]))
    if report["unexpected_metrics"]:
        print("benchmark_report unexpected_metrics=" + ",".join(report["unexpected_metrics"]))


def main() -> int:
    args = parse_args()
    policy_path = Path(args.policy).resolve()
    current_path = Path(args.current).resolve()
    policy_baseline_path, default_threshold_pct, metric_policies = load_policy(policy_path, args.profile)

    baseline_path = Path(args.baseline).resolve() if args.baseline else policy_baseline_path
    if baseline_path is None:
        raise SystemExit("baseline non definita")

    baseline = read_metrics(baseline_path)
    current = read_metrics(current_path)

    report = compare_metrics(
        baseline=baseline,
        current=current,
        metric_policies=metric_policies,
        default_threshold_pct=default_threshold_pct,
        strict_metrics=args.strict_metrics,
    )
    report.update(
        {
            "profile": args.profile,
            "policy_file": str(policy_path),
            "baseline_file": str(baseline_path),
            "current_file": str(current_path),
        }
    )

    emit_text_report(report)

    if args.report_json:
        report_path = Path(args.report_json).resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return 0 if report["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
