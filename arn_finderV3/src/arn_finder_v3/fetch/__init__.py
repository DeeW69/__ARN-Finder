"""Fetch V3 — client NCBI avec retry tenacity + fetch async aiohttp."""

from __future__ import annotations

import asyncio
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import aiohttp
from tenacity import retry, stop_after_attempt, wait_exponential, before_sleep_log

from ..config import NCBISettings

logger = logging.getLogger("arnfinderv3.fetch")

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


@dataclass
class FetchRequest:
    query: str
    db: str = "nucleotide"
    limit: int = 200
    batch_size: int = 50
    rettype: str = "fasta"
    retmode: str = "text"
    out_dir: Path = Path("data/raw")


@dataclass
class FetchResult:
    ids: list[str] = field(default_factory=list)
    fasta_path: Optional[Path] = None
    metadata_path: Optional[Path] = None
    manifest_path: Optional[Path] = None
    fetched: int = 0
    elapsed: float = 0.0


# ── Client asynchrone ─────────────────────────────────────────────────────────

class NCBIAsyncClient:
    """
    Client NCBI Entrez asynchrone (aiohttp).

    V3 vs V1/V2 :
    - aiohttp remplace requests → batch FASTA fetch en //
    - tenacity gère les retry avec backoff exponentiel
    - Rate limiter intégré via asyncio.Semaphore
    """

    def __init__(self, settings: NCBISettings) -> None:
        self._settings = settings
        # Limite le nombre de requêtes simultanées selon la clé API
        max_concurrency = 5 if settings.api_key else 2
        self._semaphore = asyncio.Semaphore(max_concurrency)

    def _base_params(self) -> dict[str, str]:
        params: dict[str, str] = {
            "email": self._settings.email,
            "tool": self._settings.tool,
        }
        if self._settings.api_key:
            params["api_key"] = self._settings.api_key
        return params

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        before_sleep=before_sleep_log(logger, logging.WARNING),
        reraise=True,
    )
    async def _get(self, session: aiohttp.ClientSession, url: str, params: dict) -> str:
        async with self._semaphore:
            async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=60)) as resp:
                resp.raise_for_status()
                return await resp.text()

    async def esearch(self, session: aiohttp.ClientSession, db: str, query: str, retmax: int) -> list[str]:
        params = {**self._base_params(), "db": db, "term": query, "retmax": str(retmax), "retmode": "json"}
        text = await self._get(session, f"{ENTREZ_BASE}/esearch.fcgi", params)
        import json
        data = json.loads(text)
        return data["esearchresult"]["idlist"]

    async def efetch_batch(self, session: aiohttp.ClientSession, ids: list[str], db: str, rettype: str, retmode: str) -> str:
        params = {**self._base_params(), "db": db, "id": ",".join(ids), "rettype": rettype, "retmode": retmode}
        return await self._get(session, f"{ENTREZ_BASE}/efetch.fcgi", params)


async def _async_fetch(req: FetchRequest, settings: NCBISettings) -> FetchResult:
    """Implémentation async interne."""
    client = NCBIAsyncClient(settings)
    result = FetchResult()
    t0 = time.perf_counter()

    req.out_dir.mkdir(parents=True, exist_ok=True)

    async with aiohttp.ClientSession() as session:
        ids = await client.esearch(session, req.db, req.query, req.limit)
        result.ids = ids
        logger.info("%d IDs récupérés", len(ids))

        # Fetch FASTA en batches concurrents
        batches = [ids[i:i + req.batch_size] for i in range(0, len(ids), req.batch_size)]
        tasks = [client.efetch_batch(session, b, req.db, req.rettype, req.retmode) for b in batches]
        responses = await asyncio.gather(*tasks)

    fasta_path = req.out_dir / "sequences.fasta"
    fasta_path.write_text("\n".join(responses), encoding="utf-8")

    result.fasta_path = fasta_path
    result.fetched = len(ids)
    result.elapsed = time.perf_counter() - t0
    logger.info("Fetch terminé en %.1fs → %s", result.elapsed, fasta_path)
    return result


def fetch_sequences(req: FetchRequest, settings: NCBISettings) -> FetchResult:
    """Point d'entrée synchrone (utilise asyncio.run en interne)."""
    return asyncio.run(_async_fetch(req, settings))
