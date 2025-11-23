"""Quick diagnostic to test NCBI calls without running the Flask server.
Run with: python test_ncbi.py
"""
from app import fetch_ncbi_genes
from Bio import Entrez


def main():
    print("Entrez.email:", Entrez.email)
    q = "hypertrophic cardiomyopathy"
    print(f"Querying NCBI for '{q}'...\n")
    results = fetch_ncbi_genes(q, retmax=10)
    print(f"Returned {len(results)} results")
    for i, r in enumerate(results[:10], start=1):
        print(f"{i}. {r}")


if __name__ == "__main__":
    main()
