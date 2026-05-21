"""Tests du client fetch V3 (mock aiohttp avec respx)."""

from __future__ import annotations

import pytest
import respx
import httpx

# NOTE: les tests réels seront implémentés avec respx pour mocker aiohttp
# Squelette posé pour forcer TDD lors du dev du module fetch


def test_placeholder_fetch_module_importable() -> None:
    """Le module fetch doit être importable sans erreur."""
    from arn_finder_v3.fetch import FetchRequest, FetchResult  # noqa: F401
    assert True
