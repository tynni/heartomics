# Heartomics
*A bioinformatics web app for exploring genomic data related to heart disease.*

---

## Overview
**Heartomics** is a Flask-based bioinformatics project that analyzes the **genomic properties of heart disease**.  
It connects to public databases (like **NCBI**) to retrieve and display genes associated with cardiovascular disease, helping visualize genetic patterns and hereditary factors.

---

## Features
- Fetches heart-disease–related gene data from **NCBI** using Biopython  
- Displays results in a clean, searchable web table  

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

### 2 Create a Virtual Environment (OPTIONAL: For MAC Users)
```bash
python3 -m venv venv
source venv/bin/activate     # Mac/Linux
venv\Scripts\activate        # Windows

### 3 Install Dependencies
```bash
pip install -r requirements.txt

### 4 Create a .env File
```bash
NCBI_EMAIL=your_email@example.com
NCBI_API_KEY=your_optional_api_key

### 5 Run the Flask App
```bash
python3 app.py
