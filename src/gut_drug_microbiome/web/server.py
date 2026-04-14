from __future__ import annotations

import json
import traceback
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler
from http.server import ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs
from urllib.parse import urlparse

from .service import DEFAULT_DEMO_RANKING_PATH
from .service import DEFAULT_INTEGRATED_PREDICTIONS_PATH
from .service import GutPredictionService


ROOT = Path(__file__).resolve().parents[3]
DEFAULT_STATIC_DIR = ROOT / "webapp" / "static"


class GutPredictionHTTPServer(ThreadingHTTPServer):
    daemon_threads = True

    def __init__(
        self,
        server_address: tuple[str, int],
        request_handler_class: type[SimpleHTTPRequestHandler],
        service: GutPredictionService,
        static_dir: str | Path,
    ) -> None:
        self.service = service
        self.static_dir = Path(static_dir)
        super().__init__(server_address, request_handler_class)


class GutPredictionRequestHandler(SimpleHTTPRequestHandler):
    server: GutPredictionHTTPServer

    def __init__(self, request, client_address, server) -> None:  # type: ignore[override]
        super().__init__(request, client_address, server, directory=str(server.static_dir))

    def _write_json(self, status: int, payload: dict[str, object]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def _write_error(self, status: int, message: str) -> None:
        self._write_json(status, {"error": message})

    def _require_query_param(self, params: dict[str, list[str]], key: str) -> str:
        values = params.get(key, [])
        if not values or not values[0].strip():
            raise ValueError(f"缺少参数: {key}")
        return values[0].strip()

    def _optional_text(self, payload: dict[str, object], key: str, default: str | None = None) -> str | None:
        if key not in payload:
            return default
        value = payload.get(key)
        if value is None:
            return default
        text = str(value).strip()
        if not text:
            return default
        if text.lower() in {"none", "null", "undefined"}:
            return default
        return text

    def do_GET(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if parsed.path.startswith("/api/"):
            try:
                query = parse_qs(parsed.query)
                if parsed.path == "/api/bootstrap":
                    self._write_json(HTTPStatus.OK, self.server.service.bootstrap())
                    return
                if parsed.path == "/api/drug-profile":
                    drug = self._require_query_param(query, "drug")
                    self._write_json(HTTPStatus.OK, self.server.service.get_drug_profile(drug))
                    return
                if parsed.path == "/api/pair-prediction":
                    drug = self._require_query_param(query, "drug")
                    microbe = self._require_query_param(query, "microbe")
                    self._write_json(HTTPStatus.OK, self.server.service.get_pair_prediction(drug, microbe))
                    return
                if parsed.path == "/api/custom-drug/pair":
                    session_id = self._require_query_param(query, "session_id")
                    microbe = self._require_query_param(query, "microbe")
                    self._write_json(
                        HTTPStatus.OK,
                        self.server.service.get_custom_pair_prediction(session_id, microbe),
                    )
                    return
                self._write_error(HTTPStatus.NOT_FOUND, f"未知接口: {parsed.path}")
            except ValueError as exc:
                self._write_error(HTTPStatus.BAD_REQUEST, str(exc))
            except Exception as exc:  # noqa: BLE001
                traceback.print_exc()
                self._write_error(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))
            return

        if parsed.path in {"", "/"}:
            self.path = "/index.html"
        else:
            self.path = parsed.path
        super().do_GET()

    def do_POST(self) -> None:  # noqa: N802
        parsed = urlparse(self.path)
        if not parsed.path.startswith("/api/"):
            self._write_error(HTTPStatus.NOT_FOUND, f"未知接口: {parsed.path}")
            return

        try:
            content_length = int(self.headers.get("Content-Length", "0"))
            raw = self.rfile.read(content_length) if content_length > 0 else b"{}"
            payload = json.loads(raw.decode("utf-8") or "{}")

            if parsed.path == "/api/custom-drug/predict":
                self._write_json(
                    HTTPStatus.OK,
                    self.server.service.predict_custom_drug(
                        drug_name=self._optional_text(payload, "drug_name"),
                        smiles=self._optional_text(payload, "smiles", default="") or "",
                        microbe_query=self._optional_text(payload, "microbe"),
                        therapeutic_class=self._optional_text(payload, "therapeutic_class"),
                        therapeutic_effect=self._optional_text(payload, "therapeutic_effect"),
                        target_species=self._optional_text(payload, "target_species", default="human") or "human",
                        human_use=bool(payload.get("human_use", True)),
                        veterinary=bool(payload.get("veterinary", False)),
                    ),
                )
                return

            if parsed.path == "/api/step3/simulate":
                self._write_json(
                    HTTPStatus.OK,
                    self.server.service.simulate_step3(
                        drug_query=self._optional_text(payload, "drug", default="") or "",
                        scenario_name=self._optional_text(payload, "scenario", default="healthy_reference")
                        or "healthy_reference",
                        community_table_path=self._optional_text(payload, "community_table_path"),
                        disease_name=self._optional_text(payload, "disease_name"),
                        n_steps=int(payload.get("n_steps", 14)),
                        initial_dose=float(payload.get("initial_dose", 1.0)),
                        repeat_dose=float(payload.get("repeat_dose", 1.0)),
                        dosing_interval=int(payload.get("dosing_interval", 1)),
                        drug_clearance_rate=float(payload.get("drug_clearance_rate", 0.12)),
                        product_clearance_rate=float(payload.get("product_clearance_rate", 0.18)),
                        metabolism_scale=float(payload.get("metabolism_scale", 0.85)),
                        effect_scale=float(payload.get("effect_scale", 0.55)),
                        ecology_strength=float(payload.get("ecology_strength", 0.20)),
                    ),
                )
                return

            if parsed.path == "/api/custom-drug/step3/simulate":
                session_id = self._optional_text(payload, "session_id", default="") or ""
                self._write_json(
                    HTTPStatus.OK,
                    self.server.service.simulate_custom_step3(
                        session_id=session_id,
                        scenario_name=self._optional_text(payload, "scenario", default="healthy_reference")
                        or "healthy_reference",
                        community_table_path=self._optional_text(payload, "community_table_path"),
                        disease_name=self._optional_text(payload, "disease_name"),
                        n_steps=int(payload.get("n_steps", 14)),
                        initial_dose=float(payload.get("initial_dose", 1.0)),
                        repeat_dose=float(payload.get("repeat_dose", 1.0)),
                        dosing_interval=int(payload.get("dosing_interval", 1)),
                        drug_clearance_rate=float(payload.get("drug_clearance_rate", 0.12)),
                        product_clearance_rate=float(payload.get("product_clearance_rate", 0.18)),
                        metabolism_scale=float(payload.get("metabolism_scale", 0.85)),
                        effect_scale=float(payload.get("effect_scale", 0.55)),
                        ecology_strength=float(payload.get("ecology_strength", 0.20)),
                    ),
                )
                return

            if parsed.path == "/api/step3/scenario-grid":
                self._write_json(
                    HTTPStatus.OK,
                    self.server.service.scenario_grid_step3(
                        drug_query=self._optional_text(payload, "drug", default="") or "",
                        community_table_path=self._optional_text(payload, "community_table_path"),
                        disease_name=self._optional_text(payload, "disease_name"),
                        n_steps=int(payload.get("n_steps", 14)),
                        initial_dose=float(payload.get("initial_dose", 1.0)),
                        repeat_dose=float(payload.get("repeat_dose", 1.0)),
                        dosing_interval=int(payload.get("dosing_interval", 1)),
                        drug_clearance_rate=float(payload.get("drug_clearance_rate", 0.12)),
                        product_clearance_rate=float(payload.get("product_clearance_rate", 0.18)),
                        metabolism_scale=float(payload.get("metabolism_scale", 0.85)),
                        effect_scale=float(payload.get("effect_scale", 0.55)),
                        ecology_strength=float(payload.get("ecology_strength", 0.20)),
                    ),
                )
                return

            if parsed.path == "/api/custom-drug/step3/scenario-grid":
                session_id = self._optional_text(payload, "session_id", default="") or ""
                self._write_json(
                    HTTPStatus.OK,
                    self.server.service.scenario_grid_custom_step3(
                        session_id=session_id,
                        community_table_path=self._optional_text(payload, "community_table_path"),
                        disease_name=self._optional_text(payload, "disease_name"),
                        n_steps=int(payload.get("n_steps", 14)),
                        initial_dose=float(payload.get("initial_dose", 1.0)),
                        repeat_dose=float(payload.get("repeat_dose", 1.0)),
                        dosing_interval=int(payload.get("dosing_interval", 1)),
                        drug_clearance_rate=float(payload.get("drug_clearance_rate", 0.12)),
                        product_clearance_rate=float(payload.get("product_clearance_rate", 0.18)),
                        metabolism_scale=float(payload.get("metabolism_scale", 0.85)),
                        effect_scale=float(payload.get("effect_scale", 0.55)),
                        ecology_strength=float(payload.get("ecology_strength", 0.20)),
                    ),
                )
                return

            self._write_error(HTTPStatus.NOT_FOUND, f"未知接口: {parsed.path}")
        except ValueError as exc:
            self._write_error(HTTPStatus.BAD_REQUEST, str(exc))
        except json.JSONDecodeError:
            self._write_error(HTTPStatus.BAD_REQUEST, "请求体不是合法 JSON。")
        except Exception as exc:  # noqa: BLE001
            traceback.print_exc()
            self._write_error(HTTPStatus.INTERNAL_SERVER_ERROR, str(exc))


def serve_web_app(
    host: str = "127.0.0.1",
    port: int = 8080,
    integrated_predictions_path: str | Path = DEFAULT_INTEGRATED_PREDICTIONS_PATH,
    demo_ranking_path: str | Path | None = DEFAULT_DEMO_RANKING_PATH,
    static_dir: str | Path = DEFAULT_STATIC_DIR,
) -> None:
    service = GutPredictionService(
        integrated_predictions_path=integrated_predictions_path,
        demo_ranking_path=demo_ranking_path,
    )
    server = GutPredictionHTTPServer(
        server_address=(host, port),
        request_handler_class=GutPredictionRequestHandler,
        service=service,
        static_dir=static_dir,
    )
    url = f"http://{host}:{port}"
    print(f"Step 1/2/3 web app is available at {url}", flush=True)
    print("Press Ctrl+C to stop the server.", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down web app server...", flush=True)
    finally:
        server.server_close()
