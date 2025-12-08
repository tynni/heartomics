# Heartomics
*A bioinformatics web app for exploring cardiovascular-genomics data from multiple databases.*

---
# Heartomics

A small Flask demo for exploring heart-related gene information across two public gene data sources.

This repository provides a simple UI to search by term or gene symbol and view results from:

- NCBI Gene (via Biopython Entrez)
- MyGene.info (public web API)

## Features

- Search NCBI Gene and MyGene by keyword or gene symbol
- View tabular results for each source (gene_id, symbol, chromosome number, and UCSC link)
- Compare overlaps (NCBI ∩ MyGene) for a query

## Requirements

- Python 3.8+
- Virtual environment recommended
- The dependencies in `requirements.txt` (install with `pip install -r requirements.txt`)

## Setup

1. Clone the repository

```bash
git clone https://github.com/yourusername/heartomics.git
cd heartomics
```

2. Create and activate a virtual environment

```bash
python3 -m venv venv
source venv/bin/activate   # macOS / Linux
```

3. Install Python dependencies

```bash
pip install -r requirements.txt
```

4. Create a `.env` file in the project root with at minimum your NCBI contact email:

```
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_optional_ncbi_api_key
SECRET_KEY=a-random-secret-for-flask-sessions
```

Notes:
- `NCBI_EMAIL` is required by NCBI Entrez utilities and should be set to an email address.
- `NCBI_API_KEY` is optional but recommended if you will do many requests.

## Running the app

Start the Flask development server:

```bash
python app.py
```

Then open your browser at `http://127.0.0.1:5000/`.

## Usage examples

- Search by gene symbol (recommended for reliable matches): `MYH7`, `TNNT2`, `APOB`, `TTN`, `MYBPC3`
- Search by keyword: `cardiomyopathy`, `heart development`, `cardiac` (may return broader results)
- Use the radio buttons to choose `NCBI`, `MyGene`, or `Compare` then submit.

### Compare behavior

- The `Compare` option computes the intersection of exact gene symbols returned by NCBI and MyGene for the given query.
- The compare view renders a table (same layout as individual results) with one row per overlapping symbol and columns:
	- `symbol` — the common HGNC symbol
	- `ncbi_gene_id`, `ncbi_name` — values from NCBI Gene
	- `mygene_id`, `mygene_name` — values from MyGene.info
	- `UCSC` — link to the UCSC Genome Browser search for the symbol (hg38)

If no exact symbol overlap is found, the app shows an informational message and the results table will be empty.

## Quick diagnostic

There is a small diagnostic helper `test_ncbi.py` that can be used to check Entrez configuration and simple queries. However, it is not necessary. It is only to test if NCBI is not connected. 

```bash
python3 app.py
python3 test_ncbi.py
```