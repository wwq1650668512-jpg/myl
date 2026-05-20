from __future__ import annotations

import argparse
import json
import math
import statistics
import sys
import threading
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gut_drug_microbiome.web.server import DEFAULT_STATIC_DIR
from gut_drug_microbiome.web.server import GutPredictionHTTPServer
from gut_drug_microbiome.web.server import GutPredictionRequestHandler
from gut_drug_microbiome.web.server import GutPredictionService
from gut_drug_microbiome.web.service import DEFAULT_INTEGRATED_PREDICTIONS_PATH


def _request_json(
    method: str,
    url: str,
    payload: dict[str, Any] | None = None,
    timeout_seconds: float = 180.0,
) -> tuple[int, float, int]:
    body: bytes | None = None
    headers = {"Accept": "application/json"}
    if payload is not None:
        body = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    request = urllib.request.Request(url=url, data=body, headers=headers, method=method.upper())
    start = time.perf_counter()
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        response_body = response.read()
        latency_ms = (time.perf_counter() - start) * 1000.0
        return response.status, latency_ms, len(response_body)


def _invoke_callable(func) -> tuple[int, float, int]:
    start = time.perf_counter()
    payload = func()
    latency_ms = (time.perf_counter() - start) * 1000.0
    body_size = len(json.dumps(payload, ensure_ascii=False))
    return 200, latency_ms, body_size


def _percentile(values: list[float], percentile: float) -> float:
    if not values:
        return math.nan
    if len(values) == 1:
        return float(values[0])
    ordered = sorted(values)
    rank = (len(ordered) - 1) * percentile
    lower = math.floor(rank)
    upper = math.ceil(rank)
    if lower == upper:
        return float(ordered[lower])
    weight = rank - lower
    return float(ordered[lower] * (1.0 - weight) + ordered[upper] * weight)


def _summarize(name: str, mode: str, latencies_ms: list[float], duration_seconds: float, failures: int) -> dict[str, Any]:
    completed = len(latencies_ms)
    requests_total = completed + failures
    throughput_rps = completed / duration_seconds if duration_seconds > 0 else math.nan
    return {
        "endpoint": name,
        "mode": mode,
        "requests_total": requests_total,
        "requests_ok": completed,
        "requests_failed": failures,
        "duration_seconds": round(duration_seconds, 4),
        "throughput_rps": round(throughput_rps, 4),
        "latency_ms_mean": round(statistics.mean(latencies_ms), 4) if latencies_ms else math.nan,
        "latency_ms_min": round(min(latencies_ms), 4) if latencies_ms else math.nan,
        "latency_ms_p50": round(_percentile(latencies_ms, 0.50), 4) if latencies_ms else math.nan,
        "latency_ms_p95": round(_percentile(latencies_ms, 0.95), 4) if latencies_ms else math.nan,
        "latency_ms_p99": round(_percentile(latencies_ms, 0.99), 4) if latencies_ms else math.nan,
        "latency_ms_max": round(max(latencies_ms), 4) if latencies_ms else math.nan,
    }


def _run_sequential(name: str, method: str, url: str, payload: dict[str, Any] | None, runs: int, timeout_seconds: float) -> dict[str, Any]:
    latencies_ms: list[float] = []
    failures = 0
    started = time.perf_counter()
    for _ in range(runs):
        try:
            status, latency_ms, _ = _request_json(method, url, payload, timeout_seconds=timeout_seconds)
            if status != 200:
                failures += 1
                continue
            latencies_ms.append(latency_ms)
        except Exception:
            failures += 1
    duration_seconds = time.perf_counter() - started
    return _summarize(name, "single_user_single_request", latencies_ms, duration_seconds, failures)


def _run_sequential_callable(name: str, func, runs: int) -> dict[str, Any]:
    latencies_ms: list[float] = []
    failures = 0
    started = time.perf_counter()
    for _ in range(runs):
        try:
            status, latency_ms, _ = _invoke_callable(func)
            if status != 200:
                failures += 1
                continue
            latencies_ms.append(latency_ms)
        except Exception:
            failures += 1
    duration_seconds = time.perf_counter() - started
    return _summarize(name, "single_user_single_request", latencies_ms, duration_seconds, failures)


def _run_concurrent(
    name: str,
    method: str,
    url: str,
    payload: dict[str, Any] | None,
    total_requests: int,
    concurrency: int,
    timeout_seconds: float,
) -> dict[str, Any]:
    latencies_ms: list[float] = []
    failures = 0
    lock = threading.Lock()

    def _worker() -> None:
        nonlocal failures
        try:
            status, latency_ms, _ = _request_json(method, url, payload, timeout_seconds=timeout_seconds)
            with lock:
                if status == 200:
                    latencies_ms.append(latency_ms)
                else:
                    failures += 1
        except Exception:
            with lock:
                failures += 1

    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = [executor.submit(_worker) for _ in range(total_requests)]
        for future in as_completed(futures):
            future.result()
    duration_seconds = time.perf_counter() - started
    return _summarize(name, f"multi_user_multi_request_c{concurrency}", latencies_ms, duration_seconds, failures)


def _run_concurrent_callable(name: str, func, total_requests: int, concurrency: int) -> dict[str, Any]:
    latencies_ms: list[float] = []
    failures = 0
    lock = threading.Lock()

    def _worker() -> None:
        nonlocal failures
        try:
            status, latency_ms, _ = _invoke_callable(func)
            with lock:
                if status == 200:
                    latencies_ms.append(latency_ms)
                else:
                    failures += 1
        except Exception:
            with lock:
                failures += 1

    started = time.perf_counter()
    with ThreadPoolExecutor(max_workers=concurrency) as executor:
        futures = [executor.submit(_worker) for _ in range(total_requests)]
        for future in as_completed(futures):
            future.result()
    duration_seconds = time.perf_counter() - started
    return _summarize(name, f"multi_user_multi_request_c{concurrency}", latencies_ms, duration_seconds, failures)


def _build_markdown(results: list[dict[str, Any]], metadata: dict[str, Any]) -> str:
    lines = [
        "# Web API 性能测试结果",
        "",
        f"- 基础地址：`{metadata['base_url']}`",
        f"- 测试药物：`{metadata['drug']}`",
        f"- 测试微生物：`{metadata['microbe']}`",
        f"- 单用户单请求轮次：`{metadata['single_runs']}`",
        f"- 多用户并发数：`{metadata['multi_concurrency']}`",
        f"- 多用户总请求数：`{metadata['multi_requests']}`",
        "",
        "| 接口 | 模式 | 成功/总数 | 吞吐 (req/s) | 平均延迟 (ms) | P50 (ms) | P95 (ms) | P99 (ms) | 最大延迟 (ms) |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in results:
        lines.append(
            "| {endpoint} | {mode} | {requests_ok}/{requests_total} | {throughput_rps} | {latency_ms_mean} | "
            "{latency_ms_p50} | {latency_ms_p95} | {latency_ms_p99} | {latency_ms_max} |".format(**row)
        )
    return "\n".join(lines) + "\n"


def _wait_until_ready(base_url: str, timeout_seconds: float = 300.0) -> None:
    started = time.perf_counter()
    while time.perf_counter() - started < timeout_seconds:
        try:
            status, _, _ = _request_json("GET", f"{base_url}/api/bootstrap", timeout_seconds=30.0)
            if status == 200:
                return
        except Exception:
            time.sleep(0.5)
    raise TimeoutError(f"Server at {base_url} was not ready within {timeout_seconds} seconds.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark local web API endpoints for single and concurrent loads.")
    parser.add_argument("--base-url", default="http://127.0.0.1:8097")
    parser.add_argument("--drug", default="Metformin hydrochloride")
    parser.add_argument("--microbe", default="Bacteroides vulgatus")
    parser.add_argument("--single-runs", type=int, default=10)
    parser.add_argument("--multi-concurrency", type=int, default=8)
    parser.add_argument("--multi-requests", type=int, default=80)
    parser.add_argument("--step3-n-steps", type=int, default=14)
    parser.add_argument("--timeout-seconds", type=float, default=180.0)
    parser.add_argument("--execution-mode", choices=["http", "service"], default="service")
    parser.add_argument("--spawn-server", action="store_true")
    parser.add_argument("--integrated-predictions", type=Path, default=None)
    parser.add_argument("--static-dir", type=Path, default=DEFAULT_STATIC_DIR)
    parser.add_argument("--output-json", type=Path, default=None)
    parser.add_argument("--output-md", type=Path, default=None)
    args = parser.parse_args()

    server: GutPredictionHTTPServer | None = None
    server_thread: threading.Thread | None = None
    service: GutPredictionService | None = None
    if args.execution_mode == "http" and args.spawn_server:
        parsed = urllib.parse.urlparse(args.base_url)
        host = parsed.hostname or "127.0.0.1"
        port = parsed.port or 8097
        service = GutPredictionService(
            integrated_predictions_path=args.integrated_predictions
            if args.integrated_predictions is not None
            else DEFAULT_INTEGRATED_PREDICTIONS_PATH
        )
        server = GutPredictionHTTPServer(
            server_address=(host, port),
            request_handler_class=GutPredictionRequestHandler,
            service=service,
            static_dir=args.static_dir,
        )
        server_thread = threading.Thread(target=server.serve_forever, daemon=True)
        server_thread.start()
        _wait_until_ready(args.base_url, timeout_seconds=300.0)
    elif args.execution_mode == "service":
        service = GutPredictionService(
            integrated_predictions_path=args.integrated_predictions
            if args.integrated_predictions is not None
            else DEFAULT_INTEGRATED_PREDICTIONS_PATH
        )

    drug_encoded = urllib.parse.quote(args.drug, safe="")
    microbe_encoded = urllib.parse.quote(args.microbe, safe="")
    endpoints = [
        {
            "name": "bootstrap",
            "method": "GET",
            "url": f"{args.base_url}/api/bootstrap",
            "payload": None,
        },
        {
            "name": "pair_prediction",
            "method": "GET",
            "url": f"{args.base_url}/api/pair-prediction?drug={drug_encoded}&microbe={microbe_encoded}",
            "payload": None,
        },
        {
            "name": "step3_simulate",
            "method": "POST",
            "url": f"{args.base_url}/api/step3/simulate",
            "payload": {
                "drug": args.drug,
                "scenario": "healthy_reference",
                "n_steps": args.step3_n_steps,
                "initial_dose": 1.0,
                "repeat_dose": 1.0,
                "dosing_interval": 1,
                "drug_clearance_rate": 0.12,
                "product_clearance_rate": 0.18,
                "metabolism_scale": 0.85,
                "effect_scale": 0.55,
                "ecology_strength": 0.20,
            },
        },
    ]

    service_endpoints = []
    if service is not None:
        service_endpoints = [
            ("bootstrap", lambda: service.bootstrap()),
            ("pair_prediction", lambda: service.get_pair_prediction(args.drug, args.microbe)),
            (
                "step3_simulate",
                lambda: service.simulate_step3(
                    drug_query=args.drug,
                    scenario_name="healthy_reference",
                    n_steps=args.step3_n_steps,
                    initial_dose=1.0,
                    repeat_dose=1.0,
                    dosing_interval=1,
                    drug_clearance_rate=0.12,
                    product_clearance_rate=0.18,
                    metabolism_scale=0.85,
                    effect_scale=0.55,
                    ecology_strength=0.20,
                ),
            ),
        ]

    results: list[dict[str, Any]] = []
    try:
        results = []
        if args.execution_mode == "http":
            for endpoint in endpoints:
                results.append(
                    _run_sequential(
                        endpoint["name"],
                        endpoint["method"],
                        endpoint["url"],
                        endpoint["payload"],
                        runs=args.single_runs,
                        timeout_seconds=args.timeout_seconds,
                    )
                )
                results.append(
                    _run_concurrent(
                        endpoint["name"],
                        endpoint["method"],
                        endpoint["url"],
                        endpoint["payload"],
                        total_requests=args.multi_requests,
                        concurrency=args.multi_concurrency,
                        timeout_seconds=args.timeout_seconds,
                    )
                )
        else:
            for endpoint_name, endpoint_callable in service_endpoints:
                results.append(_run_sequential_callable(endpoint_name, endpoint_callable, runs=args.single_runs))
                results.append(
                    _run_concurrent_callable(
                        endpoint_name,
                        endpoint_callable,
                        total_requests=args.multi_requests,
                        concurrency=args.multi_concurrency,
                    )
                )
        metadata = {
            "base_url": args.base_url,
            "drug": args.drug,
            "microbe": args.microbe,
            "single_runs": args.single_runs,
            "multi_concurrency": args.multi_concurrency,
            "multi_requests": args.multi_requests,
            "step3_n_steps": args.step3_n_steps,
            "execution_mode": args.execution_mode,
            "spawn_server": args.spawn_server,
        }
        payload = {"metadata": metadata, "results": results}

        if args.output_json is not None:
            args.output_json.parent.mkdir(parents=True, exist_ok=True)
            args.output_json.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        if args.output_md is not None:
            args.output_md.parent.mkdir(parents=True, exist_ok=True)
            args.output_md.write_text(_build_markdown(results, metadata), encoding="utf-8")

        print(json.dumps(payload, ensure_ascii=False, indent=2))
    finally:
        if server is not None:
            server.shutdown()
            server.server_close()
        if server_thread is not None:
            server_thread.join(timeout=5.0)


if __name__ == "__main__":
    main()
