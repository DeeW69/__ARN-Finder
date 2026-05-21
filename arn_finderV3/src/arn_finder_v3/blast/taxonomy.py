"""
Enrichissement taxonomique des hits BLAST — V3.2

Utilise NCBI efetch (db=taxonomy) pour obtenir la lignée complète
(kingdom / phylum / class / order / family / genus) à partir des taxid
extraits des résultats BLAST.

Les taxids sont envoyés par lots de 50 pour respecter les limites NCBI.
Fallback gracieux : si un taxid est introuvable, les champs sont laissés vides.
"""

from __future__ import annotations

import asyncio
import logging
import xml.etree.ElementTree as ET
from typing import Optional

import aiohttp
from tenacity import retry, stop_after_attempt, wait_exponential, before_sleep_log

from ..io import get_logger

logger = get_logger("arnfinderv3.blast.taxonomy")

_EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_BATCH_SIZE = 50  # max taxids par requête efetch

# Rangs taxonomiques cibles
_TARGET_RANKS = ("superkingdom", "phylum", "class", "order", "family", "genus")


@retry(
    stop=stop_after_attempt(4),
    wait=wait_exponential(multiplier=1, min=2, max=20),
    before_sleep=before_sleep_log(logger, logging.WARNING),
    reraise=True,
)
async def _fetch_taxonomy_batch(
    session: aiohttp.ClientSession,
    taxids: list[str],
    email: str,
    tool: str,
    api_key: Optional[str],
) -> dict[str, dict[str, str]]:
    """
    Récupère la lignée taxonomique pour un lot de taxids via NCBI efetch.

    Retourne : {taxid: {kingdom, phylum, class, order, family, genus}}
    """
    params: dict[str, str] = {
        "db": "taxonomy",
        "id": ",".join(taxids),
        "retmode": "xml",
        "email": email,
        "tool": tool,
    }
    if api_key:
        params["api_key"] = api_key

    async with session.get(
        _EFETCH_URL, params=params, timeout=aiohttp.ClientTimeout(total=30)
    ) as resp:
        resp.raise_for_status()
        xml_text = await resp.text()

    return _parse_taxonomy_xml(xml_text)


def _parse_taxonomy_xml(xml_text: str) -> dict[str, dict[str, str]]:
    """Parse le XML NCBI taxonomy et retourne un dict taxid → lignée."""
    result: dict[str, dict[str, str]] = {}
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError as e:
        logger.warning("Erreur parsing XML taxonomy : %s", e)
        return result

    for taxon in root.findall("Taxon"):
        taxid_el = taxon.find("TaxId")
        if taxid_el is None or not taxid_el.text:
            continue
        taxid = taxid_el.text.strip()

        lineage: dict[str, str] = {r: "" for r in _TARGET_RANKS}

        lineage_ex = taxon.find("LineageEx")
        if lineage_ex is not None:
            for node in lineage_ex.findall("Taxon"):
                rank_el = node.find("Rank")
                name_el = node.find("ScientificName")
                if rank_el is not None and name_el is not None:
                    rank = (rank_el.text or "").strip()
                    name = (name_el.text or "").strip()
                    if rank in _TARGET_RANKS:
                        lineage[rank] = name

        # Normalise "superkingdom" → "kingdom" pour les sorties
        lineage["kingdom"] = lineage.pop("superkingdom", "")
        result[taxid] = lineage

    return result


async def enrich_hits_with_taxonomy(
    hits: list[dict],
    email: str,
    tool: str = "arnfinder_v3",
    api_key: Optional[str] = None,
) -> list[dict]:
    """
    Enrichit la liste de hits BLAST avec la lignée taxonomique.

    Ajoute les champs : taxid, kingdom, phylum, class, order, family, genus.
    Les taxids manquants ou invalides donnent des champs vides.
    """
    if not hits:
        return hits

    # Collecte les taxids uniques valides
    unique_taxids = list({
        str(h.get("taxid", ""))
        for h in hits
        if h.get("taxid") and str(h["taxid"]).isdigit() and str(h["taxid"]) != "0"
    })

    if not unique_taxids:
        logger.debug("Aucun taxid valide dans les hits BLAST — enrichissement ignoré")
        _add_empty_taxonomy(hits)
        return hits

    logger.info("Enrichissement taxonomique pour %d taxids uniques", len(unique_taxids))
    lineage_map: dict[str, dict[str, str]] = {}

    # Batching + rate limiting
    async with aiohttp.ClientSession() as session:
        for i in range(0, len(unique_taxids), _BATCH_SIZE):
            batch = unique_taxids[i:i + _BATCH_SIZE]
            try:
                batch_result = await _fetch_taxonomy_batch(session, batch, email, tool, api_key)
                lineage_map.update(batch_result)
            except Exception as e:
                logger.warning("Échec efetch taxonomy (batch %d) : %s", i // _BATCH_SIZE, e)
            # Respecte le rate limit NCBI (3 req/s sans clé)
            await asyncio.sleep(0.35 if not api_key else 0.12)

    logger.info("Lignées récupérées pour %d/%d taxids", len(lineage_map), len(unique_taxids))

    # Enrichit chaque hit
    empty_lineage = {r: "" for r in ("kingdom", "phylum", "class", "order", "family", "genus")}
    for hit in hits:
        taxid = str(hit.get("taxid", ""))
        lin = lineage_map.get(taxid, empty_lineage)
        hit.update({
            "kingdom": lin.get("kingdom", ""),
            "phylum":  lin.get("phylum", ""),
            "class":   lin.get("class", ""),
            "order":   lin.get("order", ""),
            "family":  lin.get("family", ""),
            "genus":   lin.get("genus", ""),
        })

    return hits


def _add_empty_taxonomy(hits: list[dict]) -> None:
    """Ajoute des champs taxonomiques vides quand il n'y a pas de taxid."""
    for hit in hits:
        hit.update(kingdom="", phylum="", **{"class": ""}, order="", family="", genus="")
