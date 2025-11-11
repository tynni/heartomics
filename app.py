from flask import Flask, render_template, request, redirect, url_for, flash
from Bio import Entrez
import os
import pandas as pd
import requests
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

app = Flask(__name__)
app.secret_key = os.getenv("SECRET_KEY", "dev-secret")


def fetch_ncbi_genes(term: str, retmax: int = 50):
    """Query NCBI Gene for a term and return a list of dicts with gene info.

    Uses Entrez.esearch + Entrez.esummary. Returns an empty list on error.
    """
    if not Entrez.email:
        raise RuntimeError("NCBI_EMAIL not set in environment; set it in .env")

    try:
        handle = Entrez.esearch(db="gene", term=f"{term} AND Homo sapiens[Organism]", retmax=retmax)
        record = Entrez.read(handle)
        handle.close()
        ids = record.get("IdList", [])
        if not ids:
            return []

        # fetch summaries
        handle = Entrez.esummary(db="gene", id=','.join(ids))
        summaries = Entrez.read(handle)
        handle.close()

        results = []
        for s in summaries:
            results.append({
                "gene_id": s.get("Id"),
                "symbol": s.get("NomenclatureSymbol") or s.get("Name"),
                "name": s.get("Description"),
                "chromosome": s.get("Chromosome"),
            })

        return results
    except Exception as e:
        # keep failure silent for UX; caller can display message
        app.logger.exception("NCBI fetch failed")
        return []


def fetch_reactome_pathways(term: str, species: str = "Homo sapiens", size: int = 50):
    """Search Reactome ContentService for pathways matching term.

    Returns list of dicts {"stId","name","species"} or empty list on failure.
    """
    try:
        q = requests.utils.quote(term)
        url = f"https://reactome.org/ContentService/search/query?q={q}&species=Homo+sapiens&type=Pathway&size={size}"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        results = []
        for item in data.get("results", [])[:size]:
            results.append({
                "stId": item.get("stId") or item.get("dbId"),
                "name": item.get("name"),
                "species": item.get("species"),
            })
        return results
    except Exception:
        app.logger.exception("Reactome fetch failed")
        return []


def fetch_clingen_genes(term: str, size: int = 50):
    """Attempt to query ClinGen (best-effort). If the public API is not reachable,
    return an empty list.

    ClinGen APIs vary; this function attempts a couple of known endpoints and
    falls back gracefully.
    """
    try:
        # Try ClinGen Search service (best-effort).
        q = requests.utils.quote(term)
        url = f"https://search.clinicalgenome.org/api?q={q}&size={size}"
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()

        results = []
        # API shapes vary; attempt to parse expected structures
        hits = data.get("hits") or data.get("results") or []
        for h in hits[:size]:
            name = h.get("label") or h.get("symbol") or h.get("name")
            gid = h.get("id") or h.get("geneId")
            results.append({"gene_id": gid, "symbol": name, "name": h.get("description")})
        return results
    except Exception:
        app.logger.exception("ClinGen fetch failed")
        return []


@app.route("/")
def index():
    # simple index that shows a search box and tabs for each data source
    return render_template("index.html")


@app.route("/search")
def search():
    q = request.args.get("q", "heart")
    source = request.args.get("source", "ncbi")

    if source == "ncbi":
        try:
            results = fetch_ncbi_genes(q)
        except RuntimeError as re:
            flash(str(re), "danger")
            results = []
        app.logger.info(f"NCBI search for '{q}' returned {len(results)} results")
        if not results:
            flash("NCBI returned no results for that query.", "info")
        return render_template("results.html", source=source.upper(), query=q, rows=results)

    if source == "reactome":
        results = fetch_reactome_pathways(q)
        app.logger.info(f"Reactome search for '{q}' returned {len(results)} results")
        if not results:
            flash("Reactome returned no results for that query.", "info")
        return render_template("results.html", source=source.upper(), query=q, rows=results)

    if source == "clingen":
        results = fetch_clingen_genes(q)
        app.logger.info(f"ClinGen search for '{q}' returned {len(results)} results")
        if not results:
            flash("ClinGen returned no results for that query.", "info")
        return render_template("results.html", source=source.upper(), query=q, rows=results)

    if source == "compare":
        # fetch all three and compute simple overlaps by symbol/name
        ncbi = fetch_ncbi_genes(q)
        clingen = fetch_clingen_genes(q)
        reactome = fetch_reactome_pathways(q)

        app.logger.info(f"Compare for '{q}': ncbi={len(ncbi)}, clingen={len(clingen)}, reactome={len(reactome)}")

        # normalize keys to symbol/name for set operations
        ncbi_set = set([r.get("symbol") for r in ncbi if r.get("symbol")])
        clingen_set = set([r.get("symbol") for r in clingen if r.get("symbol")])
        reactome_set = set([r.get("name") for r in reactome if r.get("name")])

        overlap = {
            "ncbi_clingen": sorted(list(ncbi_set & clingen_set)),
            "ncbi_reactome": sorted(list(ncbi_set & reactome_set)),
            "clingen_reactome": sorted(list(clingen_set & reactome_set)),
        }

        return render_template("compare.html", query=q, overlap=overlap)

    # unknown source -> redirect to index
    return redirect(url_for("index"))


if __name__ == "__main__":
    # Standard/simple entrypoint: run on the default Flask host/port.
    # Keep this minimal so users can override the port externally if desired.
    app.run(debug=True)