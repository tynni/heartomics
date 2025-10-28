from flask import Flask, render_template
from Bio import Entrez
import os
import pandas as pd
from dotenv import load_dotenv

load_dotenv()

Entrez.email = os.getenv("NCBI_EMAIL")
Entrez.api_key = os.getenv("NCBI_API_KEY")

app = Flask(__name__)

if __name__ == "__main__":
    app.run(debug=True)