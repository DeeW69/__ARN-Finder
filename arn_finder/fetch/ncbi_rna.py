"""NCBI Entrez client for fetching RNA/DNA sequences."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import requests

from ..io.paths import ensure_dir, now_iso, sha256_digest, write_json, write_lines, write_text

ENTREZ_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ESEARCH_ENDPOINT = f"{ENTREZ_BASE}/esearch.fcgi"
EFETCH_ENDPOINT = f"{ENTREZ_BASE}/efetch.fcgi"
ESUMMARY_ENDPOINT = f"{ENTREZ_BASE}/esummary.fcgi"

ESEARCH_BATCH_SIZE = 500
DEFAULT_FETCH_BATCH_SIZE = 100


@dataclass(slots=True)
class FetchRequest:
    """Parameters used to query NCBI."""

    query: str
    limit: int
    out_dir: Path
    cache_dir: Path
    db: str = "nuccore"
    rettype: str = "fasta"
    retmode: str = "text"
    email: str | None = None
    tool: str = "arnfinder"
    api_key: str | None = None
    sleep: float = 0.34
    timeout: int = 30
    retries: int = 3
    batch_size: int = DEFAULT_FETCH_BATCH_SIZE


@dataclass(slots=True)
class FetchResult:
    """Summary of outputs produced by fetch."""

    ids: list[str]
    total_count: int
    fasta_path: Path
    metadata_path: Path
    ids_path: Path
    manifest_path: Path


class NCBIClient:
    """Thin wrapper around requests with caching and retries."""

    def __init__(self, request: FetchRequest, logger) -> None:
        self.request = request
        self.logger = logger
        self.session = requests.Session()

    def esearch(self, retstart: int, retmax: int) -> dict[str, Any]:
        params = {
            "db": self.request.db,
            "term": self.request.query,
            "retstart": retstart,
            "retmax": retmax,
            "retmode": "json",
            "usehistory": "y",
        }
        return self._get_json(
            ESEARCH_ENDPOINT,
            params,
            cache_kind="esearch",
            cache_suffix=".json",
            cache_key=f"{self.request.db}|{self.request.query}|{retstart}|{retmax}",
        )

    def efetch(self, ids: Sequence[str]) -> str:
        params = {
            "db": self.request.db,
            "id": ",".join(ids),
            "rettype": self.request.rettype,
            "retmode": self.request.retmode,
        }
        return self._get_text(
            EFETCH_ENDPOINT,
            params,
            cache_kind="efetch",
            cache_suffix=".fasta",
            cache_key=f"{self.request.db}|{','.join(ids)}|{self.request.rettype}|{self.request.retmode}",
        )

    def esummary(self, ids: Sequence[str]) -> dict[str, Any]:
        params = {
            "db": self.request.db,
            "id": ",".join(ids),
            "retmode": "json",
        }
        return self._get_json(
            ESUMMARY_ENDPOINT,
            params,
            cache_kind="esummary",
            cache_suffix=".json",
            cache_key=f"{self.request.db}|{','.join(ids)}",
        )

    def _request(
        self,
        url: str,
        params: dict[str, Any],
        *,
        expect_json: bool,
        cache_path: Path | None,
    ) -> str | dict[str, Any]:
        if cache_path and cache_path.exists():
            self.logger.debug("Cache hit: %s", cache_path)
            if expect_json:
                with cache_path.open("r", encoding="utf-8") as handle:
                    return json.load(handle)
            return cache_path.read_text(encoding="utf-8")

        combined_params = {**params, **self._shared_params()}
        last_exc: Exception | None = None
        for attempt in range(1, self.request.retries + 1):
            try:
                response = self.session.get(url, params=combined_params, timeout=self.request.timeout)
                response.raise_for_status()
                data = response.json() if expect_json else response.text
                if cache_path:
                    ensure_dir(cache_path.parent)
                    if expect_json:
                        write_json(cache_path, data)
                    else:
                        write_text(cache_path, data)
                if self.request.sleep > 0:
                    time.sleep(self.request.sleep)
                return data
            except requests.RequestException as exc:
                last_exc = exc
                delay = 2** (attempt - 1)
                self.logger.warning(
                    "Request to %s failed (attempt %s/%s): %s",
                    url,
                    attempt,
                    self.request.retries,
                    exc,
                )
                if attempt == self.request.retries:
                    break
                time.sleep(delay)
        assert last_exc is not None
        raise RuntimeError(f"Failed to call {url}: {last_exc}") from last_exc

    def _get_json(
        self,
        url: str,
        params: dict[str, Any],
        *,
        cache_kind: str,
        cache_suffix: str,
        cache_key: str,
    ) -> dict[str, Any]:
        cache_path = self._cache_path(cache_kind, cache_key, cache_suffix)
        return self._request(url, params, expect_json=True, cache_path=cache_path)  # type: ignore[return-value]

    def _get_text(
        self,
        url: str,
        params: dict[str, Any],
        *,
        cache_kind: str,
        cache_suffix: str,
        cache_key: str,
    ) -> str:
        cache_path = self._cache_path(cache_kind, cache_key, cache_suffix)
        return self._request(url, params, expect_json=False, cache_path=cache_path)  # type: ignore[return-value]

    def _shared_params(self) -> dict[str, Any]:
        params: dict[str, Any] = {"tool": self.request.tool}
        if self.request.email:
            params["email"] = self.request.email
        if self.request.api_key:
            params["api_key"] = self.request.api_key
        return params

    def _cache_path(self, kind: str, cache_key: str, suffix: str) -> Path:
        digest = sha256_digest(cache_key)
        path = self.request.cache_dir / kind / f"{digest}{suffix}"
        ensure_dir(path.parent)
        return path


def fetch_sequences(request: FetchRequest, logger) -> FetchResult:
    """Fetch public sequences and metadata from NCBI."""

    ensure_dir(request.out_dir)
    ensure_dir(request.cache_dir)
    client = NCBIClient(request, logger)

    ids, total_count = _collect_ids(client, request, logger)
    fasta_path = request.out_dir / "sequences.fasta"
    metadata_path = request.out_dir / "metadata.jsonl"
    ids_path = request.out_dir / "ids.txt"
    manifest_path = request.out_dir / "manifest.json"

    if not ids:
        write_text(fasta_path, "")
        write_text(metadata_path, "")
        write_lines(ids_path, [])
        manifest = _build_manifest(request, total_count, ids, fasta_path, metadata_path, ids_path)
        write_json(manifest_path, manifest)
        return FetchResult(ids=ids, total_count=total_count, fasta_path=fasta_path, metadata_path=metadata_path, ids_path=ids_path, manifest_path=manifest_path)

    logger.info("Downloading %s sequences (found %s total)", len(ids), total_count)
    _write_fasta(client, request, ids, fasta_path)
    _write_metadata(client, request, ids, metadata_path, logger)
    write_lines(ids_path, ids)

    manifest = _build_manifest(request, total_count, ids, fasta_path, metadata_path, ids_path)
    write_json(manifest_path, manifest)

    return FetchResult(
        ids=ids,
        total_count=total_count,
        fasta_path=fasta_path,
        metadata_path=metadata_path,
        ids_path=ids_path,
        manifest_path=manifest_path,
    )


def _collect_ids(client: NCBIClient, request: FetchRequest, logger) -> tuple[list[str], int]:
    ids: list[str] = []
    retstart = 0
    total_count = 0
    while len(ids) < request.limit:
        retmax = min(ESEARCH_BATCH_SIZE, request.limit - len(ids))
        response = client.esearch(retstart=retstart, retmax=retmax)
        result = response.get("esearchresult", {})
        total_count = int(result.get("count", 0))
        batch_ids: list[str] = result.get("idlist", [])
        if not batch_ids:
            break
        ids.extend(batch_ids)
        retstart += len(batch_ids)
        logger.debug("Fetched %s IDs (retstart=%s)", len(ids), retstart)
        if retstart >= total_count:
            break
    trimmed_ids = ids[: request.limit]
    logger.info("ESearch returned %s IDs (requesting %s)", len(trimmed_ids), request.limit)
    return trimmed_ids, total_count


def _write_fasta(client: NCBIClient, request: FetchRequest, ids: Sequence[str], output_path: Path) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for batch in _chunked(ids, request.batch_size):
            fasta_text = client.efetch(batch)
            if fasta_text and not fasta_text.endswith("\n"):
                fasta_text += "\n"
            handle.write(fasta_text)


def _write_metadata(
    client: NCBIClient,
    request: FetchRequest,
    ids: Sequence[str],
    output_path: Path,
    logger,
) -> None:
    with output_path.open("w", encoding="utf-8") as handle:
        for batch in _chunked(ids, request.batch_size):
            try:
                summary = client.esummary(batch)
                result = summary.get("result", {})
            except RuntimeError as exc:
                logger.warning("Metadata fetch failed for %s IDs, writing minimal entries: %s", len(batch), exc)
                result = {}
            for uid in batch:
                entry = result.get(uid, {})
                metadata = _build_metadata_row(uid, entry, request)
                handle.write(json.dumps(metadata) + "\n")


def _build_metadata_row(uid: str, entry: dict[str, Any], request: FetchRequest) -> dict[str, Any]:
    return {
        "uid": uid,
        "accession": entry.get("accessionversion") or entry.get("accession"),
        "title": entry.get("title"),
        "organism": entry.get("organism") or entry.get("taxname"),
        "slen": entry.get("slen") or entry.get("length"),
        "taxid": entry.get("taxid"),
        "source_db": request.db,
        "query": request.query,
    }


def _build_manifest(
    request: FetchRequest,
    total_count: int,
    ids: Sequence[str],
    fasta_path: Path,
    metadata_path: Path,
    ids_path: Path,
) -> dict[str, Any]:
    return {
        "query": request.query,
        "timestamp": now_iso(),
        "db": request.db,
        "limit_requested": request.limit,
        "total_found": total_count,
        "ids_retrieved": list(ids),
        "parameters": {
            "rettype": request.rettype,
            "retmode": request.retmode,
            "tool": request.tool,
            "sleep": request.sleep,
            "timeout": request.timeout,
            "retries": request.retries,
            "email_present": bool(request.email),
            "api_key_present": bool(request.api_key),
        },
        "outputs": {
            "sequences": str(fasta_path),
            "metadata": str(metadata_path),
            "ids": str(ids_path),
        },
        "cache_dir": str(request.cache_dir),
    }


def _chunked(seq: Sequence[str], size: int) -> Iterable[Sequence[str]]:
    for idx in range(0, len(seq), size):
        yield seq[idx : idx + size]
