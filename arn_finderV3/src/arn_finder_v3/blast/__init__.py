"""Client BLAST V3 — polling asynchrone avec tenacity + cache local."""

from __future__ import annotations

import asyncio
import enum
import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import aiohttp
from tenacity import retry, stop_after_attempt, wait_exponential, before_sleep_log

from ..io import ensure_dir, now_iso, write_json, write_csv, get_logger
from .parser import parse_blast_json

logger = get_logger("arnfinderv3.blast")

BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


class BlastJobStatus(enum.Enum):
    WAITING = "WAITING"
    READY = "READY"
    FAILED = "FAILED"


class BlastClientError(RuntimeError):
    pass


@dataclass(slots=True)
class BlastRequest:
    input_fasta: Path
    out_dir: Path
    program: str = "blastn"
    db: str = "nt"
    max_hits: int = 10
    poll_interval: int = 10
    timeout: int = 600
    email: str = ""
    tool: str = "arnfinder_v3"
    api_key: Optional[str] = None
    enrich_taxonomy: bool = True      # V3.2 — enrichissement taxonomique des hits


@dataclass
class BlastResult:
    hits_jsonl: Path
    summary_csv: Path
    manifest_path: Path
    total_queries: int
    total_hits: int
    elapsed: float = 0.0


# ── Client asynchrone ─────────────────────────────────────────────────────────

class BlastAsyncClient:
    """Client BLAST asynchrone avec retry tenacity sur submit/poll/fetch."""

    def __init__(self, email: str, tool: str, api_key: Optional[str] = None) -> None:
        if not email:
            raise BlastClientError("Email NCBI requis pour BLAST.")
        self._email = email
        self._tool = tool
        self._api_key = api_key

    def _base_params(self) -> dict:
        p = {"EMAIL": self._email, "TOOL": self._tool}
        if self._api_key:
            p["API_KEY"] = self._api_key
        return p

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        before_sleep=before_sleep_log(logger, logging.WARNING),
        reraise=True,
    )
    async def _post(self, session: aiohttp.ClientSession, data: dict) -> str:
        async with session.post(BASE_URL, data=data, timeout=aiohttp.ClientTimeout(total=60)) as resp:
            resp.raise_for_status()
            return await resp.text()

    @retry(
        stop=stop_after_attempt(5),
        wait=wait_exponential(multiplier=1, min=2, max=30),
        before_sleep=before_sleep_log(logger, logging.WARNING),
        reraise=True,
    )
    async def _get(self, session: aiohttp.ClientSession, params: dict) -> str:
        async with session.get(BASE_URL, params=params, timeout=aiohttp.ClientTimeout(total=60)) as resp:
            resp.raise_for_status()
            return await resp.text()

    async def submit(self, session: aiohttp.ClientSession, sequence: str, program: str, db: str) -> tuple[str, Optional[int]]:
        payload = {**self._base_params(), "CMD": "Put", "PROGRAM": program, "DATABASE": db, "QUERY": sequence}
        text = await self._post(session, payload)
        rid: Optional[str] = None
        rtoe: Optional[int] = None
        for line in text.splitlines():
            if "RID" in line and "=" in line:
                rid = line.partition("=")[2].strip()
            elif "RTOE" in line and "=" in line:
                try:
                    rtoe = int(line.partition("=")[2].strip())
                except ValueError:
                    pass
        if not rid:
            raise BlastClientError("Pas de RID dans la réponse BLAST.")
        return rid, rtoe

    async def poll(self, session: aiohttp.ClientSession, rid: str) -> BlastJobStatus:
        params = {**self._base_params(), "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"}
        text = await self._get(session, params)
        for line in text.splitlines():
            if "=" in line:
                key, _, val = line.partition("=")
                if key.strip().upper() == "STATUS":
                    normalized = val.strip().upper()
                    if normalized == "READY":
                        return BlastJobStatus.READY
                    if normalized == "FAILED":
                        return BlastJobStatus.FAILED
        return BlastJobStatus.WAITING

    async def fetch_results(self, session: aiohttp.ClientSession, rid: str) -> str:
        params = {**self._base_params(), "CMD": "Get", "RID": rid, "FORMAT_TYPE": "JSON2"}
        return await self._get(session, params)

    async def run_query(
        self, session: aiohttp.ClientSession,
        sequence: str, record_id: str, program: str, db: str,
        poll_interval: int, timeout: int, max_hits: int,
    ) -> list[dict]:
        """Submit → poll → fetch pour une seule séquence."""
        rid, rtoe = await self.submit(session, sequence, program, db)
        initial_sleep = max(rtoe or poll_interval, poll_interval)
        logger.debug("RID=%s RTOE=%s pour %s", rid, rtoe, record_id)
        await asyncio.sleep(initial_sleep)

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            status = await self.poll(session, rid)
            if status == BlastJobStatus.READY:
                break
            if status == BlastJobStatus.FAILED:
                logger.warning("BLAST échoué pour %s (RID=%s)", record_id, rid)
                return []
            await asyncio.sleep(poll_interval)
        else:
            logger.warning("Timeout BLAST pour %s (RID=%s)", record_id, rid)
            return []

        raw = await self.fetch_results(session, rid)
        hits = parse_blast_json(raw, max_hits, len(sequence))
        for h in hits:
            h["query_id"] = record_id
        return hits


# ── Orchestrateur ─────────────────────────────────────────────────────────────

async def _async_blast(request: BlastRequest) -> BlastResult:
    from ..io import read_fasta_records

    out_dir = ensure_dir(request.out_dir)
    hits_path = out_dir / "blast_hits.jsonl"
    summary_path = out_dir / "blast_summary.csv"
    manifest_path = out_dir / "blast_manifest.json"
    raw_dir = ensure_dir(out_dir / "raw")

    t0 = time.perf_counter()
    client = BlastAsyncClient(request.email, request.tool, request.api_key)
    all_hits: list[dict] = []
    total_queries = 0

    # Séquentiel pour respecter les limites NCBI BLAST (pas de // sur BLAST)
    async with aiohttp.ClientSession() as session:
        for header, sequence in read_fasta_records(request.input_fasta):
            record_id = header.split(" ", 1)[0]
            total_queries += 1
            logger.info("BLAST [%s/%s] %s", total_queries, "?", record_id)
            hits = await client.run_query(
                session, sequence, record_id,
                request.program, request.db,
                request.poll_interval, request.timeout, request.max_hits,
            )
            all_hits.extend(hits)
            (raw_dir / f"{record_id}.json").write_text(
                json.dumps({"record_id": record_id, "hits": hits}, indent=2), encoding="utf-8"
            )

    # V3.2 — Enrichissement taxonomique des hits BLAST
    if request.enrich_taxonomy and all_hits:
        try:
            from .taxonomy import enrich_hits_with_taxonomy
            all_hits = await enrich_hits_with_taxonomy(
                all_hits, request.email, request.tool, request.api_key
            )
            logger.info("Enrichissement taxonomique terminé pour %d hits", len(all_hits))
        except Exception as e:
            logger.warning("Enrichissement taxonomique échoué (hits non enrichis) : %s", e)
            # Ajoute des champs vides pour garder le schéma cohérent
            for h in all_hits:
                h.setdefault("taxid", "")
                for col in ("kingdom", "phylum", "class", "order", "family", "genus"):
                    h.setdefault(col, "")
    else:
        for h in all_hits:
            h.setdefault("taxid", "")
            for col in ("kingdom", "phylum", "class", "order", "family", "genus"):
                h.setdefault(col, "")

    with hits_path.open("w", encoding="utf-8") as f:
        for hit in all_hits:
            f.write(json.dumps(hit) + "\n")

    summary_rows = [
        [
            h["query_id"], h["accession"], h["organism"],
            h.get("taxid", ""),
            f"{h['pident']:.3f}", f"{h['evalue']}", f"{h['bitscore']:.1f}", f"{h['coverage']:.3f}",
            h.get("kingdom", ""), h.get("phylum", ""), h.get("class", ""),
            h.get("order", ""), h.get("family", ""), h.get("genus", ""),
        ]
        for h in all_hits
    ]
    write_csv(
        summary_path,
        ["query_id", "accession", "organism", "taxid",
         "pident", "evalue", "bitscore", "coverage",
         "kingdom", "phylum", "class", "order", "family", "genus"],
        summary_rows,
    )

    elapsed = time.perf_counter() - t0
    manifest = {
        "timestamp": now_iso(),
        "params": {
            "program": request.program, "db": request.db,
            "max_hits": request.max_hits, "poll_interval": request.poll_interval, "timeout": request.timeout,
        },
        "inputs": {"fasta": str(request.input_fasta)},
        "total_queries": total_queries,
        "total_hits": len(all_hits),
        "elapsed_s": round(elapsed, 1),
        "outputs": {"hits_jsonl": str(hits_path), "summary_csv": str(summary_path)},
    }
    write_json(manifest_path, manifest)
    logger.info("BLAST terminé : %s queries, %s hits en %.1fs", total_queries, len(all_hits), elapsed)

    return BlastResult(
        hits_jsonl=hits_path,
        summary_csv=summary_path,
        manifest_path=manifest_path,
        total_queries=total_queries,
        total_hits=len(all_hits),
        elapsed=elapsed,
    )


def run_blast(request: BlastRequest) -> BlastResult:
    """Point d'entrée synchrone."""
    return asyncio.run(_async_blast(request))
