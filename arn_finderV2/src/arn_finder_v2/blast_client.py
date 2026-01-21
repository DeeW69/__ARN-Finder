"""Thin wrapper around the NCBI BLAST HTTP API."""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass
from typing import Optional

import requests

logger = logging.getLogger(__name__)


class BlastClientError(RuntimeError):
    """Base error raised for BLAST client failures."""


@dataclass
class BlastStatus:
    state: str
    message: Optional[str] = None


class BlastClient:
    """Minimal BLAST HTTP client with throttling."""

    BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    def __init__(
        self,
        email: str,
        tool: str,
        api_key: Optional[str] = None,
        rate_limit_seconds: float = 1.0,
        session: Optional[requests.Session] = None,
    ) -> None:
        if not email:
            raise ValueError("email is required")
        if not tool:
            raise ValueError("tool is required")
        self.email = email
        self.tool = tool
        self.api_key = api_key
        self.rate_limit_seconds = max(rate_limit_seconds, 0.0)
        self.session = session or requests.Session()
        self._last_request_at = 0.0

    def submit(self, query_seq: str, program: str, db: str) -> str:
        """Submit a query sequence and return the RID."""
        logger.info("Submitting sequence (%d bp) to BLAST %s/%s", len(query_seq), program, db)
        data = {
            "CMD": "Put",
            "PROGRAM": program,
            "DATABASE": db,
            "QUERY": query_seq,
            "EMAIL": self.email,
            "TOOL": self.tool,
        }
        if self.api_key:
            data["API_KEY"] = self.api_key

        response_text = self._post(self.BASE_URL, data=data)
        rid = self._extract_rid(response_text)
        if not rid:
            raise BlastClientError("Impossible de récupérer le RID depuis la réponse BLAST.")
        return rid

    def poll(self, rid: str, poll_seconds: int, timeout_seconds: int) -> BlastStatus:
        """Poll SearchInfo until completion or timeout."""
        deadline = time.monotonic() + timeout_seconds
        while True:
            status = self._check_status(rid)
            if status.state in {"READY", "FAILED"}:
                return status
            if time.monotonic() >= deadline:
                return BlastStatus("TIMEOUT", "Polling timeout reached.")
            logger.info("RID %s still processing, waiting %ss", rid, poll_seconds)
            time.sleep(max(poll_seconds, 1))

    def fetch(self, rid: str, format_type: str) -> str:
        """Retrieve the BLAST result in the requested format."""
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": format_type,
        }
        if self.api_key:
            params["API_KEY"] = self.api_key
        logger.info("Fetching BLAST results for %s (%s)", rid, format_type)
        return self._get(self.BASE_URL, params=params)

    # Internal helpers -------------------------------------------------

    def _check_status(self, rid: str) -> BlastStatus:
        params = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo",
        }
        if self.api_key:
            params["API_KEY"] = self.api_key

        text = self._get(self.BASE_URL, params=params)
        state = None
        message = None
        for line in text.splitlines():
            line = line.strip()
            if line.startswith("Status="):
                state = line.split("=", 1)[1].strip().upper()
            elif line.startswith("Message="):
                message = line.split("=", 1)[1].strip()

        if not state:
            raise BlastClientError(f"Réponse SearchInfo inattendue pour RID {rid}: {text}")

        return BlastStatus(state, message)

    def _extract_rid(self, response_text: str) -> Optional[str]:
        for line in response_text.splitlines():
            if "RID" in line:
                parts = line.split("=", 1)
                if len(parts) == 2:
                    return parts[1].strip()
        return None

    def _post(self, url: str, data: dict) -> str:
        self._respect_rate_limit()
        response = self.session.post(url, data=data, timeout=120)
        response.raise_for_status()
        self._last_request_at = time.monotonic()
        return response.text

    def _get(self, url: str, params: dict) -> str:
        self._respect_rate_limit()
        response = self.session.get(url, params=params, timeout=120)
        response.raise_for_status()
        self._last_request_at = time.monotonic()
        return response.text

    def _respect_rate_limit(self) -> None:
        now = time.monotonic()
        elapsed = now - self._last_request_at
        remaining = self.rate_limit_seconds - elapsed
        if remaining > 0:
            time.sleep(remaining)
