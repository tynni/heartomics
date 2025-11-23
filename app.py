from flask import Flask, render_template, request, redirect, url_for, flash
from Bio import Entrez
import os
import requests
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

app = Flask(__name__)
app.secret_key = os.getenv("SECRET_KEY", "dev-secret")

# Fetch gene information from NCBI Gene database
def fetch_ncbi_genes(term: str, retmax: int = 50):
    if not Entrez.email:
        raise RuntimeError("NCBI_EMAIL not set in environment; set it in .env")

    try:
        qraw = term.strip()

        if qraw and qraw.isalnum() and qraw.upper() == qraw and len(qraw) <= 8:
            search_term = f"{qraw}[sym] AND Homo sapiens[Organism]"
            esearch_retmax = max(retmax, 50)
        else:
            search_term = f"{term} AND Homo sapiens[Organism]"
            esearch_retmax = max(retmax, 100)

        handle = Entrez.esearch(db="gene", term=search_term, retmax=esearch_retmax)
        record = Entrez.read(handle)
        handle.close()

        ids = record.get("IdList", [])
        if not ids:
            return []

        handle = Entrez.esummary(db="gene", id=",".join(ids))
        summaries = Entrez.read(handle)
        handle.close()

        docs = summaries if isinstance(summaries, list) else summaries.get(
            "DocumentSummarySet", {}
        ).get("DocumentSummary", [])

        if not isinstance(docs, list):
            docs = [docs]

        results = []
        for i, s in enumerate(docs):
            try:
                symbol = s.get("NomenclatureSymbol") or s.get("Name")
            except Exception:
                symbol = None

            try:
                name = s.get("Description") or s.get("Title") or s.get("Summary")
            except Exception:
                name = None

            try:
                chrom = s.get("Chromosome")
            except Exception:
                chrom = None

            gid = None
            if ids and i < len(ids):
                gid = ids[i]

            ucsc_url = None
            if symbol:
                ucsc_url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={requests.utils.quote(symbol)}"

            results.append({
                "gene_id": gid or (s.get("Id") if hasattr(s, "get") else None),
                "symbol": symbol,
                "name": name,
                "chromosome": chrom,
                "ucsc_url": ucsc_url,
            })

        return results

    except Exception:
        app.logger.exception("NCBI fetch failed")
        return []


# Fetch gene information from MyGene.info
def fetch_mygene(term: str, size: int = 50):
    try:
        url = (
            "https://mygene.info/v3/query"
            f"?q={requests.utils.quote(term)}"
            "&species=human"
            f"&size={size}"
        )

        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()

        hits = data.get("hits", [])

        results = []
        for h in hits[:size]:
            symbol = h.get("symbol")
            ucsc_url = None
            if symbol:
                ucsc_url = f"https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position={requests.utils.quote(symbol)}"

            results.append({
                "gene_id": h.get("_id"),
                "symbol": symbol,
                "name": h.get("name"),
                "chromosome": h.get("chrom"),
                "ucsc_url": ucsc_url,
            })
        return results

    except Exception:
        app.logger.exception("MyGene fetch failed")
        return []

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/search")
def search():
    q = request.args.get("q", "heart")
    source = request.args.get("source", "ncbi")

    if source == "ncbi":
        results = fetch_ncbi_genes(q)
        if not results:
            flash("NCBI returned no results.", "info")
        return render_template("results.html", source="NCBI", query=q, rows=results)

    if source == "mygene":
        results = fetch_mygene(q)
        if not results:
            flash("MyGene returned no results.", "info")
        return render_template("results.html", source="MyGene", query=q, rows=results)

    if source == "compare":
        ncbi = fetch_ncbi_genes(q)
        mg = fetch_mygene(q)

        ncbi_map = {r.get("symbol"): r for r in ncbi if r.get("symbol")}
        mg_map = {r.get("symbol"): r for r in mg if r.get("symbol")}

        symbols = sorted(set(ncbi_map.keys()) & set(mg_map.keys()))

        rows = []
        for sym in symbols:
            n = ncbi_map.get(sym, {})
            m = mg_map.get(sym, {})
            rows.append({
                "symbol": sym,
                "ncbi_gene_id": n.get("gene_id"),
                "ncbi_name": n.get("name"),
                "mygene_id": m.get("gene_id"),
                "mygene_name": m.get("name"),
                "ucsc_url": n.get("ucsc_url") or m.get("ucsc_url"),
            })

        if not rows:
            flash("No overlapping symbols found between NCBI and MyGene.", "info")

        return render_template("results.html", source="Compare (NCBI ∩ MyGene)", query=q, rows=rows)

    return redirect(url_for("index"))

if __name__ == "__main__":
    app.run(debug=True)