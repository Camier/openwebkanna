from __future__ import annotations

import re
from typing import Any, cast

try:
    from rdkit import Chem
except ImportError:  # pragma: no cover - optional runtime dependency
    Chem = None

from services.multimodal_retrieval_api.adapter import (
    build_one_collection_backend_info,
    load_qdrant_client,
    load_one_collection_router,
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
    probe_exact_chemistry_runtime,
)
from scripts.rag.local_evidence_store import load_local_evidence_store


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
            "sampled_page_points_missing_figure_records": 0,
            "sampled_page_point_ids_missing_figure_records": [],
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

        missing_ids: list[str] = []
        for point in points:
            payload = getattr(point, "payload", None) or {}
            figure_records = payload.get("figure_records")
            if isinstance(figure_records, list):
                continue
            point_id = getattr(point, "id", None)
            missing_ids.append(str(point_id) if point_id is not None else "")

        report["sampled_page_points_missing_figure_records"] = len(missing_ids)
        report["sampled_page_point_ids_missing_figure_records"] = missing_ids

        if missing_ids:
            report["status"] = "degraded"
            report["warnings"].append(
                "Sampled page points in rag_evidence are missing figure_records; native evidence coverage is incomplete."
            )
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
        exact_chemistry_warnings = probe_exact_chemistry_runtime()

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
            "exact_chemistry": self._lane_status_report(
                enabled=True,
                warnings=exact_chemistry_warnings,
                metadata={
                    "collection_name": self.service_settings.rag_collection_name,
                    "object_type": "molecule",
                    "match_fields": ["inchikey", "canonical_smiles"],
                    "review_status": "parsed",
                },
            ),
        }

        warnings = [
            *config_warnings,
            *text_warnings,
            *visual_warnings,
            *exact_chemistry_warnings,
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
                "exact_chemistry": self._lane_capability_report(
                    configured=True,
                    warnings=exact_chemistry_warnings,
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
        router_result = load_one_collection_router().retrieve(
            query=request.query,
            top_k=request.top_k,
        )
        backend = build_one_collection_backend_info(
            service_settings=self.service_settings,
        )
        raw_text_result = router_result.raw_text_result or {}
        candidate_hits = self._build_candidate_hits(
            exact_hits=router_result.exact_hits,
            raw_text_result=raw_text_result,
        )
        reranked_hits = self._build_reranked_hits(router_result=router_result)
        evidence_objects, evidence_diagnostics = (
            self._build_evidence_objects_from_router_result(
                router_result=router_result,
                query=request.query,
            )
        )
        diagnostics = dict(raw_text_result.get("diagnostics") or {})
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
        diagnostics["query_type"] = router_result.query_type
        diagnostics["normalized_query"] = router_result.normalized_query
        diagnostics["lanes_used"] = router_result.lanes_used
        diagnostics["exact_chemistry_collection"] = (
            self.service_settings.rag_collection_name
        )
        diagnostics["exact_chemistry_hit_count"] = len(router_result.exact_hits)
        diagnostics["visual_lane_collection"] = (
            self.service_settings.rag_collection_name
        )
        diagnostics["visual_lane_vector_name"] = (
            self.service_settings.vision_vector_name
        )
        diagnostics["visual_hit_count"] = len(router_result.visual_hits)
        diagnostics["fused_hit_count"] = len(router_result.fused_hits)
        diagnostics.update(evidence_diagnostics)
        if router_result.warnings:
            diagnostics.setdefault("warnings", []).extend(router_result.warnings)
        response_mode = str(raw_text_result.get("mode") or request.mode or "hybrid")
        return RetrieveResponse(
            query=raw_text_result.get("query") or request.query,
            mode=response_mode,
            top_k=request.top_k,
            backend=backend,
            candidate_hits=candidate_hits,
            reranked_hits=reranked_hits,
            evidence_objects=evidence_objects,
            diagnostics=diagnostics,
        )

    def _normalize_hit(self, payload: dict[str, Any]) -> RetrievalEvidenceHit:
        return RetrievalEvidenceHit(**payload)

    def _build_candidate_hits(
        self,
        *,
        exact_hits: list[dict[str, Any]],
        raw_text_result: dict[str, Any],
    ) -> list[RetrievalEvidenceHit]:
        hits: list[RetrievalEvidenceHit] = []
        seen: set[str] = set()

        for hit in exact_hits:
            normalized = self._router_hit_as_retrieval_hit(hit)
            if normalized.point_id in seen:
                continue
            seen.add(normalized.point_id)
            hits.append(normalized)

        for hit in raw_text_result.get("candidate_hits") or []:
            normalized = self._normalize_hit(hit)
            if normalized.point_id in seen:
                continue
            seen.add(normalized.point_id)
            hits.append(normalized)
        return hits

    def _build_reranked_hits(self, *, router_result: Any) -> list[RetrievalEvidenceHit]:
        return [
            self._router_hit_as_retrieval_hit(hit) for hit in router_result.fused_hits
        ]

    def _router_hit_as_retrieval_hit(self, hit: dict[str, Any]) -> RetrievalEvidenceHit:
        payload = hit.get("payload", {}) if isinstance(hit.get("payload"), dict) else {}
        return self._normalize_hit(
            {
                "point_id": str(
                    hit.get("point_id")
                    or hit.get("id")
                    or payload.get("point_id")
                    or ""
                ),
                "rank": int(hit.get("rank") or 1),
                "score": self._payload_float(hit.get("score")) or 0.0,
                "stage": str(
                    hit.get("match_type") or payload.get("match_type") or "router_fused"
                ),
                "payload": payload,
                "doc_id": hit.get("document_id")
                or self._payload_string(payload, "document_id")
                or self._payload_string(payload, "doc_id"),
                "title": self._payload_string(payload, "document_title")
                or self._payload_string(payload, "title"),
                "page_number": hit.get("page_number")
                or self._payload_int(payload, "page_number"),
                "source_pdf_sha256": self._payload_string(payload, "source_pdf_sha256"),
                "citation_key": self._payload_string(payload, "citation_key"),
            }
        )

    def _build_evidence_objects_from_router_result(
        self,
        *,
        router_result: Any,
        query: str,
    ) -> tuple[list[RetrievalEvidenceObject], dict[str, Any]]:
        query_smiles = self._query_smiles(query)
        evidence_objects: list[RetrievalEvidenceObject] = []
        seen: set[str] = set()
        page_hits_without_native_figures = 0
        page_hit_ids_without_native_figures: list[str] = []
        page_hits_with_local_fallback_figures = 0

        for exact_hit in router_result.exact_hits:
            molecule_evidence = self._molecule_evidence_from_exact_hit(exact_hit)
            if molecule_evidence.evidence_id not in seen:
                seen.add(molecule_evidence.evidence_id)
                evidence_objects.append(molecule_evidence)
            for parent_evidence in self._expand_parent_evidence_from_exact_hit(
                exact_hit
            ):
                if parent_evidence.evidence_id in seen:
                    continue
                seen.add(parent_evidence.evidence_id)
                evidence_objects.append(parent_evidence)

        for hit in router_result.fused_hits:
            payload = (
                hit.get("payload", {}) if isinstance(hit.get("payload"), dict) else {}
            )
            object_type = str(
                hit.get("object_type") or payload.get("object_type") or "page"
            )
            if object_type in {"figure", "chemical_block"}:
                figure_evidence = self._figure_evidence_from_visual_hit(hit)
                if figure_evidence.evidence_id in seen:
                    continue
                seen.add(figure_evidence.evidence_id)
                evidence_objects.append(figure_evidence)
                continue

            if object_type == "molecule":
                molecule_evidence = self._molecule_evidence_from_exact_hit(hit)
                if molecule_evidence.evidence_id not in seen:
                    seen.add(molecule_evidence.evidence_id)
                    evidence_objects.append(molecule_evidence)
                for parent_evidence in self._expand_parent_evidence_from_exact_hit(hit):
                    if parent_evidence.evidence_id in seen:
                        continue
                    seen.add(parent_evidence.evidence_id)
                    evidence_objects.append(parent_evidence)
                continue

            page_hit = self._router_hit_as_retrieval_hit(hit)
            page_evidence = self._page_evidence_from_hit(page_hit)
            if page_evidence.evidence_id not in seen:
                seen.add(page_evidence.evidence_id)
                evidence_objects.append(page_evidence)

            native_figures = self._select_qdrant_figures(
                hit=page_hit,
                query_smiles=query_smiles,
            )
            if native_figures:
                for figure_evidence in native_figures:
                    if figure_evidence.evidence_id in seen:
                        continue
                    seen.add(figure_evidence.evidence_id)
                    evidence_objects.append(figure_evidence)
                continue
            page_hits_without_native_figures += 1
            page_hit_ids_without_native_figures.append(page_hit.point_id)
            local_figures = self._select_local_figures(
                hit=page_hit,
                query_smiles=query_smiles,
            )
            if local_figures:
                page_hits_with_local_fallback_figures += 1
                for figure_evidence in local_figures:
                    if figure_evidence.evidence_id in seen:
                        continue
                    seen.add(figure_evidence.evidence_id)
                    evidence_objects.append(figure_evidence)

        diagnostics = {
            "page_hits_without_native_figures": page_hits_without_native_figures,
            "page_hit_ids_without_native_figures": page_hit_ids_without_native_figures,
            "page_hits_with_local_fallback_figures": page_hits_with_local_fallback_figures,
        }
        if page_hits_without_native_figures:
            diagnostics["native_figure_payload_warning"] = (
                "Some page hits in rag_evidence did not include figure_records. "
                "Local extraction fallback was used when matching evidence was available."
            )
        return evidence_objects, diagnostics

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
        payload = hit.get("payload", {}) if isinstance(hit.get("payload"), dict) else {}
        point_id = str(hit.get("id") or payload.get("point_id") or "")
        figure_id = self._payload_string(payload, "figure_id") or point_id
        object_type = self._payload_string(payload, "object_type")
        evidence_type = (
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
            smiles=self._payload_list_strings(payload, "smiles"),
            canonical_smiles=self._payload_list_strings(payload, "canonical_smiles"),
            smiles_review_status=self._payload_list_strings(
                payload, "smiles_review_status"
            ),
            smiles_backends=self._payload_list_strings(payload, "smiles_backends"),
            smiles_confidences=self._payload_list_numbers(
                payload, "smiles_confidences"
            ),
            figure_kind=self._payload_string(payload, "figure_kind"),
        )

    def _molecule_evidence_from_exact_hit(
        self, hit: dict[str, Any]
    ) -> RetrievalEvidenceObject:
        payload = hit.get("payload", {}) if isinstance(hit.get("payload"), dict) else {}
        point_id = str(hit.get("id") or payload.get("point_id") or "")
        molecule_id = self._payload_string(payload, "molecule_id") or point_id
        page_number = self._payload_int(payload, "page_number")
        page_index = self._payload_int(payload, "page_index")
        if page_index is None:
            page_index = self._page_index_from_number(page_number)

        raw_smiles = self._payload_string(payload, "raw_smiles")
        canonical_smiles = self._payload_string(payload, "canonical_smiles")
        review_status = self._payload_string(payload, "review_status")
        source_text = self._payload_string(payload, "source_text")
        formula = self._payload_string(payload, "formula")
        backend = self._payload_string(payload, "backend")
        confidence = self._payload_float(payload.get("confidence"))

        text_fragments = [
            value for value in [source_text, formula, canonical_smiles] if value
        ]
        return RetrievalEvidenceObject(
            evidence_id=f"molecule:{molecule_id}",
            evidence_type="molecule",
            source="qdrant",
            doc_id=self._payload_string(payload, "document_id") or point_id,
            title=self._payload_string(payload, "document_title")
            or self._payload_string(payload, "title"),
            page_number=page_number,
            page_index=page_index,
            point_id=point_id,
            parent_point_id=self._payload_string(payload, "parent_id"),
            score=self._payload_float(hit.get("score")),
            figure_id=self._payload_string(payload, "figure_id"),
            molecule_id=molecule_id,
            block_id=self._payload_string(payload, "block_id"),
            block_type=self._payload_string(payload, "block_type"),
            text=" | ".join(text_fragments) if text_fragments else None,
            image_path=self._payload_string(payload, "image_uri")
            or self._payload_string(payload, "image_path"),
            bbox=self._payload_list_numbers(payload, "bbox"),
            smiles=[raw_smiles] if raw_smiles else [],
            canonical_smiles=[canonical_smiles] if canonical_smiles else [],
            smiles_review_status=[review_status] if review_status else [],
            smiles_backends=[backend] if backend else [],
            smiles_confidences=[confidence] if confidence is not None else [],
            figure_kind=self._payload_string(payload, "figure_kind"),
        )

    def _expand_parent_evidence_from_exact_hit(
        self, hit: dict[str, Any]
    ) -> list[RetrievalEvidenceObject]:
        payload = hit.get("payload", {}) if isinstance(hit.get("payload"), dict) else {}
        point_id = str(hit.get("id") or payload.get("point_id") or "")
        doc_id = self._payload_string(payload, "document_id") or point_id
        title = self._payload_string(payload, "document_title") or self._payload_string(
            payload, "title"
        )
        page_number = self._payload_int(payload, "page_number")
        page_index = self._payload_int(payload, "page_index")
        if page_index is None:
            page_index = self._page_index_from_number(page_number)

        linked: list[RetrievalEvidenceObject] = []

        figure_id = self._payload_string(payload, "figure_id")
        figure_text = self._payload_string(
            payload, "caption_text"
        ) or self._payload_string(payload, "source_text")
        image_path = self._payload_string(payload, "image_uri") or self._payload_string(
            payload, "image_path"
        )
        if figure_id or image_path or figure_text:
            linked.append(
                RetrievalEvidenceObject(
                    evidence_id=f"figure:{figure_id or point_id}",
                    evidence_type="figure",
                    source="qdrant",
                    doc_id=doc_id,
                    title=title,
                    page_number=page_number,
                    page_index=page_index,
                    point_id=figure_id or point_id,
                    parent_point_id=self._payload_string(payload, "page_id"),
                    score=max(self._payload_float(hit.get("score")) or 0.0, 0.98),
                    figure_id=figure_id,
                    caption_text=self._payload_string(payload, "caption_text"),
                    text=figure_text,
                    image_path=image_path,
                    bbox=self._payload_list_numbers(payload, "bbox"),
                    figure_kind=self._payload_string(payload, "figure_kind"),
                )
            )

        page_id = self._payload_string(payload, "page_id")
        page_text = self._payload_string(payload, "page_text") or self._payload_string(
            payload, "chunk_text"
        )
        if page_id or page_number is not None:
            linked.append(
                RetrievalEvidenceObject(
                    evidence_id=f"page:{page_id or doc_id}:{page_number or 0}",
                    evidence_type="page",
                    source="qdrant",
                    doc_id=doc_id,
                    title=title,
                    page_number=page_number,
                    page_index=page_index,
                    point_id=page_id,
                    score=max(self._payload_float(hit.get("score")) or 0.0, 0.97),
                    text=page_text,
                )
            )

        return linked

    def _select_qdrant_figures(
        self,
        *,
        hit: RetrievalEvidenceHit,
        query_smiles: str | None,
    ) -> list[RetrievalEvidenceObject]:
        figure_records = self._payload_list(hit.payload, "figure_records")
        if not figure_records:
            return []

        figures: list[RetrievalEvidenceObject] = []
        for index, record in enumerate(figure_records):
            if not isinstance(record, dict):
                continue
            figure_smiles_records = self._payload_list(record, "smiles_records")
            smiles = []
            canonical_smiles = []
            smiles_backends = []
            smiles_confidences = []
            smiles_review_status = []
            for smiles_record in figure_smiles_records:
                if not isinstance(smiles_record, dict):
                    continue
                status = smiles_record.get("status")
                if not isinstance(status, str) or status.strip() != "parsed":
                    continue
                smiles_value = smiles_record.get("smiles")
                if isinstance(smiles_value, str) and smiles_value.strip():
                    smiles.append(smiles_value.strip())
                canonical_value = smiles_record.get("canonical_smiles")
                if isinstance(canonical_value, str) and canonical_value.strip():
                    canonical_smiles.append(canonical_value.strip())
                backend_value = smiles_record.get("backend")
                if isinstance(backend_value, str) and backend_value.strip():
                    smiles_backends.append(backend_value.strip())
                confidence = smiles_record.get("confidence")
                parsed_confidence = self._payload_float(confidence)
                if parsed_confidence is not None:
                    smiles_confidences.append(parsed_confidence)
                smiles_review_status.append("parsed")

            if query_smiles is not None and not smiles and not canonical_smiles:
                continue

            if query_smiles is not None and not self._figure_smiles_match_records(
                query_smiles=query_smiles,
                smiles_records=smiles,
                canonical_smiles_records=canonical_smiles,
            ):
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
                    smiles=smiles,
                    canonical_smiles=canonical_smiles,
                    smiles_backends=smiles_backends,
                    smiles_confidences=smiles_confidences,
                    smiles_review_status=smiles_review_status,
                    figure_kind=self._payload_string(record, "figure_kind"),
                )
            )
        return figures

    def _select_local_figures(
        self,
        *,
        hit: RetrievalEvidenceHit,
        query_smiles: str | None,
    ) -> list[RetrievalEvidenceObject]:
        document = self._resolve_local_document(hit)
        if document is None:
            return []

        figures: list[RetrievalEvidenceObject] = []
        for figure in document.get_page_figures(self._page_index_from_hit(hit)):
            smiles = []
            canonical_smiles = []
            smiles_backends = []
            smiles_confidences = []
            smiles_review_status = []

            for smiles_record in figure.smiles_records:
                if smiles_record.status != "parsed":
                    continue
                if smiles_record.smiles and smiles_record.smiles.strip():
                    smiles.append(smiles_record.smiles.strip())
                if (
                    smiles_record.canonical_smiles
                    and smiles_record.canonical_smiles.strip()
                ):
                    canonical_smiles.append(smiles_record.canonical_smiles.strip())
                if smiles_record.backend and smiles_record.backend.strip():
                    smiles_backends.append(smiles_record.backend.strip())
                if smiles_record.confidence is not None:
                    smiles_confidences.append(float(smiles_record.confidence))
                smiles_review_status.append("parsed")

            if query_smiles is not None:
                if not smiles and not canonical_smiles:
                    continue
                if not self._figure_smiles_match_records(
                    query_smiles=query_smiles,
                    smiles_records=smiles,
                    canonical_smiles_records=canonical_smiles,
                ):
                    continue

            figures.append(
                RetrievalEvidenceObject(
                    evidence_id=f"{figure.evidence_type}:{figure.evidence_id}",
                    evidence_type=figure.evidence_type,
                    source="local_extraction",
                    doc_id=document.doc_id,
                    title=figure.title or hit.title,
                    page_number=hit.page_number,
                    page_index=figure.page_index,
                    point_id=hit.point_id,
                    parent_point_id=hit.point_id,
                    rank=hit.rank,
                    score=hit.score,
                    source_pdf_sha256=hit.source_pdf_sha256,
                    citation_key=hit.citation_key,
                    figure_id=figure.figure_id,
                    block_id=figure.block_id,
                    block_type=figure.block_type,
                    caption_text=figure.caption_text,
                    text=figure.text or figure.caption_text,
                    image_asset_id=figure.image_asset_id,
                    image_path=figure.image_path,
                    bbox=list(figure.bbox),
                    smiles=smiles,
                    canonical_smiles=canonical_smiles,
                    smiles_review_status=smiles_review_status,
                    smiles_backends=smiles_backends,
                    smiles_confidences=smiles_confidences,
                )
            )
        return figures

    def _resolve_local_document(self, hit: RetrievalEvidenceHit) -> Any | None:
        local_store = load_local_evidence_store()
        candidate_ids = [
            hit.doc_id,
            hit.payload.get("document_id") if isinstance(hit.payload, dict) else None,
            hit.payload.get("doc_id") if isinstance(hit.payload, dict) else None,
            hit.source_pdf_sha256,
        ]
        for candidate_id in candidate_ids:
            if isinstance(candidate_id, str):
                document = local_store.get_document(candidate_id)
                if document is not None:
                    return document

        candidate_titles = [
            hit.title,
            hit.payload.get("document_title")
            if isinstance(hit.payload, dict)
            else None,
            hit.payload.get("title") if isinstance(hit.payload, dict) else None,
        ]
        for candidate_title in candidate_titles:
            if not isinstance(candidate_title, str):
                continue
            matches = local_store.get_documents_by_normalized_title(candidate_title)
            if matches:
                return matches[0]
        return None

    def _page_index_from_hit(self, hit: RetrievalEvidenceHit) -> int | None:
        payload_page_index = (
            self._payload_int(hit.payload, "page_index")
            if isinstance(hit.payload, dict)
            else None
        )
        if payload_page_index is not None:
            return payload_page_index
        return self._page_index_from_number(hit.page_number)

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

    def _canonicalize_smiles(self, value: str | None) -> str | None:
        if Chem is None or not isinstance(value, str):
            return None
        candidate = value.strip()
        if (
            not candidate
            or len(candidate) > 512
            or any(char.isspace() for char in candidate)
        ):
            return None
        molecule = Chem.MolFromSmiles(candidate)
        if molecule is None:
            return None
        return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)

    def _query_smiles(self, value: str | None) -> str | None:
        return self._canonicalize_smiles(value)

    def _figure_smiles_match_records(
        self,
        *,
        query_smiles: str,
        smiles_records: list[str],
        canonical_smiles_records: list[str],
    ) -> bool:
        candidates: set[str] = set()
        for value in canonical_smiles_records:
            canonical = self._canonicalize_smiles(value)
            if canonical is not None:
                candidates.add(canonical)
        for value in smiles_records:
            canonical = self._canonicalize_smiles(value)
            if canonical is not None:
                candidates.add(canonical)
        return query_smiles in candidates

    def _normalize_text(self, value: str | None) -> str | None:
        if not isinstance(value, str):
            return None
        normalized = value.lower().replace("—", "-").replace("–", "-")
        normalized = re.sub(r"[^a-z0-9]+", " ", normalized)
        normalized = " ".join(normalized.split())
        return normalized or None

    def _normalize_doi(self, value: str | None) -> str | None:
        if not isinstance(value, str):
            return None
        normalized = value.strip().lower()
        normalized = normalized.removeprefix("https://doi.org/")
        normalized = normalized.removeprefix("http://doi.org/")
        normalized = normalized.removeprefix("doi:")
        return normalized or None

    def _summarize_evidence_types(
        self,
        evidence_objects: list[RetrievalEvidenceObject],
    ) -> dict[str, int]:
        counts: dict[str, int] = {}
        for item in evidence_objects:
            counts[item.evidence_type] = counts.get(item.evidence_type, 0) + 1
        return counts
