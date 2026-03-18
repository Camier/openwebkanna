from __future__ import annotations

from pydantic import ValidationError
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import JSONResponse
from starlette.routing import Route

from services.multimodal_retrieval_api.contracts import (
    RetrieveRequest,
)
from services.multimodal_retrieval_api.service import MultimodalRetrievalService


service = MultimodalRetrievalService()


async def health(_: Request) -> JSONResponse:
    return JSONResponse(service.health())


async def ready(_: Request) -> JSONResponse:
    payload = service.ready().model_dump()
    status_code = 200 if payload["status"] == "ready" else 503
    return JSONResponse(payload, status_code=status_code)


async def retrieve(request: Request) -> JSONResponse:
    try:
        payload = await request.json()
        parsed = RetrieveRequest.model_validate(payload)
    except ValidationError as exc:
        return JSONResponse(
            {"detail": exc.errors()},
            status_code=422,
        )
    except Exception as exc:
        return JSONResponse({"detail": str(exc)}, status_code=400)

    try:
        response = service.retrieve(parsed)
        return JSONResponse(response.model_dump())
    except Exception as exc:
        return JSONResponse({"detail": str(exc)}, status_code=500)


app = Starlette(
    routes=[
        Route("/health", health, methods=["GET"]),
        Route("/ready", ready, methods=["GET"]),
        Route("/api/v1/retrieve", retrieve, methods=["POST"]),
    ]
)
