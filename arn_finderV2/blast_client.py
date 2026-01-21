"""Thin wrapper around the NCBI BLAST Common URL API."""

from __future__ import annotations

import enum
import logging
from typing import Any, Dict, Tuple

import requests

LOGGER = logging.getLogger(__name__)
BASE_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


class BlastClientError(RuntimeError):
    """Raised when the BLAST API returns an unexpected payload."""


class BlastJobStatus(enum.Enum):
    """Possible job statuses returned by FORMAT_OBJECT=SearchInfo."""

    WAITING = "WAITING"
    READY = "READY"
    FAILED = "FAILED"


def parse_search_info(text: str) -> Dict[str, str]:
    """Parse a SearchInfo payload into a dictionary."""

    info: Dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        info[key.strip()] = value.strip()
    return info


def submit(sequence: str, program: str, db: str, **params: Any) -> Tuple[str, int | None]:
    """Submit a BLAST job and return (RID, RTOE)."""

    session: requests.Session | None = params.pop("session", None)
    email = params.pop("email", None)
    tool = params.pop("tool", None)
    api_key = params.pop("api_key", None)
    timeout = params.pop("timeout", 60)
    if not email:
        raise BlastClientError("NCBI email contact is required for BLAST submissions.")
    if not tool:
        raise BlastClientError("NCBI tool identifier is required for BLAST submissions.")
    payload: Dict[str, Any] = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": db,
        "QUERY": sequence,
        "EMAIL": email,
        "TOOL": tool,
    }
    optional = {k.upper(): v for k, v in params.items()}
    payload.update(optional)
    if api_key:
        payload["API_KEY"] = api_key
    sess = session or requests.Session()
    response = sess.post(BASE_URL, data=payload, timeout=timeout)
    response.raise_for_status()
    rid = None
    rtoe = None
    for line in response.text.splitlines():
        if "RID" in line:
            _, _, rid_value = line.partition("=")
            rid = rid_value.strip()
        elif "RTOE" in line:
            _, _, rtoe_value = line.partition("=")
            try:
                rtoe = int(rtoe_value.strip())
            except ValueError:
                LOGGER.debug("Unable to parse RTOE from line: %s", line)
    if not rid:
        raise BlastClientError("BLAST submission did not return a RID.")
    return rid, rtoe


def check(rid: str, **params: Any) -> BlastJobStatus:
    """Poll the SearchInfo endpoint for job status."""

    session: requests.Session | None = params.pop("session", None)
    api_key = params.pop("api_key", None)
    timeout = params.pop("timeout", 30)
    payload = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_OBJECT": "SearchInfo",
    }
    payload.update({k.upper(): v for k, v in params.items()})
    if api_key:
        payload["API_KEY"] = api_key
    sess = session or requests.Session()
    response = sess.get(BASE_URL, params=payload, timeout=timeout)
    response.raise_for_status()
    info = parse_search_info(response.text)
    status_str = info.get("Status") or info.get("STATUS")
    if not status_str:
        LOGGER.debug("SearchInfo payload missing Status: %s", response.text)
        return BlastJobStatus.WAITING
    normalized = status_str.strip().upper()
    if normalized == BlastJobStatus.READY.value:
        return BlastJobStatus.READY
    if normalized == BlastJobStatus.FAILED.value:
        return BlastJobStatus.FAILED
    return BlastJobStatus.WAITING


def fetch(rid: str, format_type: str = "JSON2", **params: Any) -> str:
    """Fetch BLAST results in the requested format."""

    session: requests.Session | None = params.pop("session", None)
    api_key = params.pop("api_key", None)
    timeout = params.pop("timeout", 60)
    payload = {"CMD": "Get", "RID": rid, "FORMAT_TYPE": format_type}
    payload.update({k.upper(): v for k, v in params.items()})
    if api_key:
        payload["API_KEY"] = api_key
    sess = session or requests.Session()
    response = sess.get(BASE_URL, params=payload, timeout=timeout)
    response.raise_for_status()
    return response.text


__all__ = ["BlastClientError", "BlastJobStatus", "submit", "check", "fetch", "parse_search_info"]
