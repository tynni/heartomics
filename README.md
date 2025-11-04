# Heartomics
*A bioinformatics web app for exploring cardiovascular-genomics data from multiple databases.*

---

## Overview
**Heartomics** is a Flask-based bioinformatics tool that retrieves and displays heart-disease–related genetic information from **three public biological databases**:

- **NCBI Gene (Entrez)** – heart-related genes
- **ClinGen** – curated gene-disease associations
- **Reactome** – cardiovascular pathways

Users can select a database via tabs/buttons, view results in a clean table, and explore a **comparison view** to see overlap across sources.

---

## Features
- Pulls heart-disease–associated genetic data from:
  - **NCBI Gene (via Biopython)**
  - **ClinGen API**
  - **Reactome Pathway API**
- Separate views/tabs for each data source
- Combined comparison view
- Searchable, structured tables
---

## Tech Stack
- **Python 3**
- **Flask**
- **Biopython**
- **pandas**
- **python-dotenv**

---

## Setup & Installation

### 1 Clone the repository
```bash
git clone https://github.com/yourusername/heartomics.git
cd heartomics
```
### 2 Create a Virtual Environment (OPTIONAL: For MAC Users)
```bash
python3 -m venv venv
source venv/bin/activate     # Mac/Linux
venv\Scripts\activate        # Windows
```
### 3 Install Dependencies
```bash
pip install -r requirements.txt
```
### 4 Create a .env File
```bash
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_optional_api_key
```
### 5 Run the Flask App
```bash
python3 app.py
```