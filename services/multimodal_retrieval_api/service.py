from __future__ import annotations

from typing import Any, cast

from services.multimodal_retrieval_api.adapter import (
    build_one_collection_backend_info,
    load_qdrant_client,
    load_retrieval_lanes,
)
from services.multimodal_retrieval_api.contracts import (
    EvidenceType,
    RetrievalEvidenceObject,
    ReadyResponse,
    RetrieveRequest,
    RetrieveResponse,
    RetrievalEvidenceHit,
)
from services.multimodal_retrieval_api.settings import ServiceSettings
from services.rag import (
    ColModernVBERTConfig,
    probe_colmodernvbert_runtime,
)


class MultimodalRetrievalService:
    def __init__(self) -> None:
        self.service_settings = ServiceSettings.load()

    def health(self) -> dict[str, str]:
        return {
            "status": "ok",
            "service": self.service_settings.app_name,
        }

    def ready(self) -> ReadyResponse:
        warnings: list[str] = []
        try:
            qdrant_client = load_qdrant_client()
            collection_info = qdrant_client.get_collection(
                self.service_settings.rag_collection_name
            )
            completeness = self._collection_completeness_report(qdrant_client)
            qdrant_status = {
                "ok": True,
                "url": self.service_settings.qdrant_url,
                "collection_name": self.service_settings.rag_collection_name,
                "points_count": getattr(collection_info, "points_count", None),
                "collection_completeness": completeness,
            }
            backend = build_one_collection_backend_info(
                service_settings=self.service_settings,
            )
            runtime = self._runtime_readiness_report()
            warnings.extend(runtime["warnings"])
            warnings.extend(completeness["warnings"])
            status = "ready" if not warnings else "degraded"
            return ReadyResponse(
                status=status,
                service=self.service_settings.app_name,
                backend=backend,
                qdrant=qdrant_status,
                runtime=runtime,
                warnings=warnings,
            )
        except Exception as exc:
            warnings.append(str(exc))
            return ReadyResponse(
                status="degraded",
                service=self.service_settings.app_name,
                backend=None,
                qdrant=None,
                runtime=None,
                warnings=warnings,
            )

    def _collection_completeness_report(self, qdrant_client: Any) -> dict[str, Any]:
        report: dict[str, Any] = {
            "status": "ready",
            "sample_limit": 100,
            "sampled_page_points": 0,
            "warnings": [],
        }
        try:
            from qdrant_client.http import models
        except ImportError:
            report["status"] = "degraded"
            report["warnings"].append(
                "Qdrant filter models are unavailable, so collection completeness could not be probed."
            )
            return report

        points, _ = qdrant_client.scroll(
            collection_name=self.service_settings.rag_collection_name,
            scroll_filter=models.Filter(
                must=[
                    models.FieldCondition(
                        key="object_type",
                        match=models.MatchValue(value="page"),
                    )
                ]
            ),
            with_payload=True,
            with_vectors=False,
            limit=report["sample_limit"],
        )
        report["sampled_page_points"] = len(points)
        return report

    def _runtime_readiness_report(self) -> dict[str, Any]:
        text_warnings: list[str] = []
        if not self.service_settings.text_query_model_path:
            text_warnings.append(
                "Dense text query model is not configured. "
                "Set MULTIMODAL_RETRIEVAL_API_TEXT_QUERY_MODEL_PATH or provide "
                "NEMOTRON_MODEL_PATH via MULTIMODAL_RETRIEVAL_API_COMPAT_ENV_FILE."
            )

        config_warnings: list[str] = []
        compat_env_file = self.service_settings.compat_env_file
        if compat_env_file is not None and not compat_env_file.exists():
            config_warnings.append(
                f"Configured compatibility env file does not exist: {compat_env_file}"
            )
        visual_warnings = probe_colmodernvbert_runtime(
            ColModernVBERTConfig().model_name,
        )

        lanes = {
            "text": self._lane_status_report(
                enabled=True,
                warnings=text_warnings,
                metadata={
                    "collection_name": self.service_settings.rag_collection_name,
                    "dense_vector_name": self.service_settings.text_dense_vector_name,
                    "sparse_vector_name": self.service_settings.text_sparse_vector_name,
                    "dense_model_path": self.service_settings.text_query_model_path,
                    "sparse_model_name": self.service_settings.text_sparse_model_name,
                    "device": self.service_settings.text_query_device,
                    "attn_implementation": self.service_settings.text_query_attn_implementation,
                },
            ),
            "visual": self._lane_status_report(
                enabled=True,
                warnings=visual_warnings,
                metadata={
                    "collection_name": self.service_settings.rag_collection_name,
                    "vector_name": self.service_settings.vision_vector_name,
                    "model_name": ColModernVBERTConfig().model_name,
                },
            ),
        }

        warnings = [
            *config_warnings,
            *text_warnings,
            *visual_warnings,
        ]
        return {
            "status": "ready" if not warnings else "degraded",
            "config": {
                "status": "ready" if not config_warnings else "degraded",
                "warnings": config_warnings,
            },
            "capabilities": {
                "text": self._lane_capability_report(
                    configured=bool(self.service_settings.text_query_model_path),
                    warnings=text_warnings,
                ),
                "visual": self._lane_capability_report(
                    configured=bool(ColModernVBERTConfig().model_name),
                    warnings=visual_warnings,
                ),
            },
            "lanes": lanes,
            "warnings": warnings,
        }

    def _lane_status_report(
        self,
        *,
        enabled: bool,
        warnings: list[str],
        metadata: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        return {
            "enabled": enabled,
            "status": "ready" if not warnings else "degraded",
            "warnings": list(warnings),
            "metadata": metadata or {},
        }

    def _lane_capability_report(
        self,
        *,
        configured: bool,
        warnings: list[str],
    ) -> dict[str, Any]:
        return {
            "configured": configured,
            "usable": configured and not warnings,
            "warning_count": len(warnings),
        }

    def retrieve(self, request: RetrieveRequest) -> RetrieveResponse:
        lanes = load_retrieval_lanes()
        text_lane = lanes["text"]
        visual_lane = lanes.get("visual")

        text_result = text_lane(query=request.query, top_k=request.top_k)

        raw_text_diagnostics = text_result.get("diagnostics") or {}

        candidate_hits = self._build_candidate_hits(text_result)
        reranked_hits = self._build_reranked_hits(text_result)

        visual_hits: list[dict[str, Any]] = []
        if visual_lane is not None:
            try:
                visual_result = visual_lane(query=request.query, top_k=request.top_k)
                visual_hits = visual_result if isinstance(visual_result, list) else []
            except Exception:
                visual_hits = []

        evidence_objects, evidence_diagnostics = self._build_evidence_objects(
            text_result=text_result,
            visual_hits=visual_hits,
        )

        diagnostics: dict[str, Any] = dict(raw_text_diagnostics)
        diagnostics["collection_report"] = {
            "collection_name": self.service_settings.rag_collection_name,
            "dense_vector_name": self.service_settings.text_dense_vector_name,
            "sparse_vector_name": self.service_settings.text_sparse_vector_name,
            "vision_vector_name": self.service_settings.vision_vector_name,
        }
        diagnostics["evidence_object_count"] = len(evidence_objects)
        diagnostics["evidence_object_types"] = self._summarize_evidence_types(
            evidence_objects
        )
        diagnostics["text_lane_hit_count"] = len(text_result.get("reranked_hits") or [])
        diagnostics["visual_lane_hit_count"] = len(visual_hits)
        diagnostics.update(evidence_diagnostics)

        backend = build_one_collection_backend_info(
            service_settings=self.service_settings,
        )
        response_mode = str(text_result.get("mode") or request.mode or "hybrid")
        return RetrieveResponse(
            query=request.query,
            mode=response_mode,
            top_k=request.top_k,
            backend=backend,
            candidate_hits=candidate_hits,
            reranked_hits=reranked_hits,
            evidence_objects=evidence_objects,
            diagnostics=diagnostics,
        )

    def _build_candidate_hits(
        self,
        text_result: dict[str, Any],
    ) -> list[RetrievalEvidenceHit]:
        hits: list[RetrievalEvidenceHit] = []
        seen: set[str] = set()
        for hit in text_result.get("candidate_hits") or []:
            normalized = self._normalize_hit(hit)
            if normalized.point_id in seen:
                continue
            seen.add(normalized.point_id)
            hits.append(normalized)
        return hits

    def _build_reranked_hits(
        self,
        text_result: dict[str, Any],
    ) -> list[RetrievalEvidenceHit]:
        hits: list[RetrievalEvidenceHit] = []
        seen: set[str] = set()
        for hit in text_result.get("reranked_hits") or []:
            normalized = self._normalize_hit(hit)
            if normalized.point_id in seen:
                continue
            seen.add(normalized.point_id)
            hits.append(normalized)
        return hits

    def _normalize_hit(self, hit: dict[str, Any]) -> RetrievalEvidenceHit:
        payload = hit.get("payload") or {}
        point_id = (
            hit.get("point_id")
            or hit.get("id")
            or payload.get("point_id")
            or payload.get("page_id")
            or ""
        )
        return RetrievalEvidenceHit(
            point_id=str(point_id),
            rank=int(hit.get("rank") or 1),
            score=float(hit.get("score") or 0.0),
            stage=str(hit.get("stage") or hit.get("match_type") or "text_hybrid"),
            doc_id=hit.get("doc_id")
            or self._payload_string(payload, "document_id")
            or self._payload_string(payload, "doc_id"),
            title=hit.get("title")
            or self._payload_string(payload, "document_title")
            or self._payload_string(payload, "title"),
            page_number=hit.get("page_number")
            or self._payload_int(payload, "page_number"),
            source_pdf_sha256=self._payload_string(payload, "source_pdf_sha256"),
            citation_key=self._payload_string(payload, "citation_key"),
            payload=payload if isinstance(payload, dict) else {},
        )

    def _build_evidence_objects(
        self,
        *,
        text_result: dict[str, Any],
        visual_hits: list[dict[str, Any]],
    ) -> tuple[list[RetrievalEvidenceObject], dict[str, Any]]:
        evidence_objects: list[RetrievalEvidenceObject] = []
        seen: set[str] = set()

        for hit in text_result.get("reranked_hits") or []:
            page_hit = self._normalize_hit(hit)
            page_evidence = self._page_evidence_from_hit(page_hit)
            if page_evidence.evidence_id not in seen:
                seen.add(page_evidence.evidence_id)
                evidence_objects.append(page_evidence)

            native_figures = self._select_qdrant_figures(hit=page_hit)
            for figure_evidence in native_figures:
                if figure_evidence.evidence_id not in seen:
                    seen.add(figure_evidence.evidence_id)
                    evidence_objects.append(figure_evidence)

        for hit in visual_hits:
            payload = hit.get("payload") or {}
            object_type = str(
                hit.get("object_type") or payload.get("object_type") or "figure"
            )
            if object_type in {"figure", "chemical_block"}:
                figure_evidence = self._figure_evidence_from_visual_hit(hit)
                if figure_evidence.evidence_id not in seen:
                    seen.add(figure_evidence.evidence_id)
                    evidence_objects.append(figure_evidence)

        return evidence_objects, {}

    def _page_evidence_from_hit(
        self, hit: RetrievalEvidenceHit
    ) -> RetrievalEvidenceObject:
        text = self._extract_page_text(hit.payload)
        return RetrievalEvidenceObject(
            evidence_id=f"page:{hit.point_id}",
            evidence_type="page",
            source="qdrant",
            doc_id=hit.doc_id or hit.source_pdf_sha256 or hit.point_id,
            title=hit.title,
            page_number=hit.page_number,
            page_index=self._page_index_from_number(hit.page_number),
            point_id=hit.point_id,
            rank=hit.rank,
            score=hit.score,
            source_pdf_sha256=hit.source_pdf_sha256,
            citation_key=hit.citation_key,
            text=text,
        )

    def _figure_evidence_from_visual_hit(
        self, hit: dict[str, Any]
    ) -> RetrievalEvidenceObject:
        payload = hit.get("payload") or {}
        point_id = str(hit.get("id") or payload.get("point_id") or "")
        figure_id = self._payload_string(payload, "figure_id") or point_id
        object_type = self._payload_string(payload, "object_type")
        evidence_type: EvidenceType = (
            "chemical_block" if object_type == "chemical_block" else "figure"
        )
        return RetrievalEvidenceObject(
            evidence_id=f"{evidence_type}:{figure_id}",
            evidence_type=evidence_type,
            source="qdrant",
            doc_id=self._payload_string(payload, "document_id") or point_id,
            title=self._payload_string(payload, "document_title")
            or self._payload_string(payload, "title"),
            page_number=self._payload_int(payload, "page_number"),
            page_index=self._payload_int(payload, "page_index")
            or self._page_index_from_number(self._payload_int(payload, "page_number")),
            point_id=point_id,
            parent_point_id=self._payload_string(payload, "parent_id"),
            score=self._payload_float(hit.get("score")),
            figure_id=self._payload_string(payload, "figure_id"),
            block_id=self._payload_string(payload, "block_id"),
            block_type=self._payload_string(payload, "block_type"),
            caption_text=self._payload_string(payload, "caption_text"),
            text=(
                self._payload_string(payload, "text")
                or self._payload_string(payload, "caption_text")
                or self._payload_string(payload, "ocr_text")
            ),
            image_path=self._payload_string(payload, "image_uri")
            or self._payload_string(payload, "image_path"),
            bbox=self._payload_list_numbers(payload, "bbox"),
            figure_kind=self._payload_string(payload, "figure_kind"),
        )

    def _select_qdrant_figures(
        self,
        *,
        hit: RetrievalEvidenceHit,
    ) -> list[RetrievalEvidenceObject]:
        figure_records = self._payload_list(hit.payload, "figure_records")
        if not figure_records:
            return []

        figures: list[RetrievalEvidenceObject] = []
        for index, record in enumerate(figure_records):
            if not isinstance(record, dict):
                continue

            evidence_type = record.get("evidence_type")
            block_id = (
                self._payload_string(record, "block_id") or f"{hit.point_id}:{index}"
            )
            if not isinstance(evidence_type, str) or evidence_type not in {
                "page",
                "figure",
                "chemical_block",
            }:
                evidence_type = (
                    "chemical_block" if "/ChemicalBlock/" in block_id else "figure"
                )

            page_number = self._payload_int(record, "page_number")
            if page_number is None:
                page_number = hit.page_number
            page_index = self._payload_int(record, "page_index")
            if page_index is None:
                page_index = self._page_index_from_number(page_number)

            evidence_id = record.get("figure_id") or block_id
            if not isinstance(evidence_id, str) or not evidence_id.strip():
                evidence_id = f"{evidence_type}:{block_id}"
            typed_evidence_type = cast(EvidenceType, evidence_type)

            figures.append(
                RetrievalEvidenceObject(
                    evidence_id=f"{evidence_type}:{evidence_id}",
                    evidence_type=typed_evidence_type,
                    source="qdrant",
                    doc_id=hit.doc_id
                    or hit.source_pdf_sha256
                    or self._payload_string(record, "document_id")
                    or hit.point_id,
                    title=self._payload_string(record, "title") or hit.title,
                    page_number=page_number,
                    page_index=page_index,
                    point_id=hit.point_id,
                    parent_point_id=hit.point_id,
                    rank=hit.rank,
                    score=hit.score,
                    source_pdf_sha256=hit.source_pdf_sha256,
                    citation_key=hit.citation_key,
                    figure_id=self._payload_string(record, "figure_id"),
                    block_id=block_id,
                    block_type=self._payload_string(record, "block_type"),
                    caption_text=self._payload_string(record, "caption_text"),
                    text=self._payload_string(record, "text")
                    or self._payload_string(record, "caption_text"),
                    image_asset_id=self._payload_string(record, "image_asset_id"),
                    image_path=self._payload_string(record, "image_uri")
                    or self._payload_string(record, "image_path"),
                    bbox=self._payload_list_numbers(record, "bbox"),
                    figure_kind=self._payload_string(record, "figure_kind"),
                )
            )
        return figures

    def _extract_page_text(self, payload: dict[str, Any]) -> str | None:
        for key in ("chunk_text", "text", "content", "page_text"):
            value = payload.get(key)
            if isinstance(value, str) and value.strip():
                return value
        return None

    def _page_index_from_number(self, page_number: int | None) -> int | None:
        if page_number is None:
            return None
        if page_number <= 0:
            return page_number
        return page_number - 1

    def _payload_string(self, payload: dict[str, Any], key: str) -> str | None:
        value = payload.get(key)
        return value if isinstance(value, str) and value.strip() else None

    def _payload_int(self, payload: dict[str, Any], key: str) -> int | None:
        value = payload.get(key)
        if value is None or value == "":
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _payload_list(self, payload: dict[str, Any], key: str) -> list[Any]:
        value = payload.get(key)
        return value if isinstance(value, list) else []

    def _payload_list_numbers(self, payload: dict[str, Any], key: str) -> list[float]:
        numbers: list[float] = []
        for value in self._payload_list(payload, key):
            parsed = self._payload_float(value)
            if parsed is not None:
                numbers.append(parsed)
        return numbers

    def _payload_list_strings(self, payload: dict[str, Any], key: str) -> list[str]:
        values: list[str] = []
        for value in self._payload_list(payload, key):
            if isinstance(value, str) and value.strip():
                values.append(value.strip())
        return values

    def _payload_float(self, value: Any) -> float | None:
        if value is None or value == "":
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _summarize_evidence_types(
        self,
        evidence_objects: list[RetrievalEvidenceObject],
    ) -> dict[str, int]:
        counts: dict[str, int] = {}
        for item in evidence_objects:
            counts[item.evidence_type] = counts.get(item.evidence_type, 0) + 1
        return counts
