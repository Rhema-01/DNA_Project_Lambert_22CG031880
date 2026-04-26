# DNA & RNA Sequence Analyzer
**CSC 442 — Project 2**  
Department of Computer Science | Faculty of Physical Sciences | 2024/2025

---

## Project Structure

```
dna_project/
├── app.py              ← Flask web application (routes, UniProt API)
├── biology.py          ← Core biology engine (detection, transcription, translation)
├── requirements.txt    ← Python dependencies
├── Procfile            ← For Render/Railway hosting
├── templates/
│   └── index.html      ← Single-page web UI
└── static/
    ├── css/style.css   ← Stylesheet
    └── js/main.js      ← Frontend logic (drag/drop, rendering, API calls)
```

---

## How to Run (VS Code)

1. Open the `dna_project` folder in VS Code
2. Open the terminal (`Ctrl + `` `)
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
4. Run the app:
   ```bash
   python app.py
   ```
5. Open your browser to: **http://127.0.0.1:5000**

---

## Features

| Feature | How it works |
|---|---|
| Sequence input | Type/paste in text area, upload file, or drag & drop |
| Sequence detection | Automatically identifies DNA, RNA, or invalid sequences with plain-English explanation |
| DNA strand type | Choose non-template or template strand (affects transcription) |
| Transcription | Produces mRNA from DNA; explains each step |
| Translation | Reads mRNA codons, maps each to amino acid with a full codon table |
| Amino acids | Displays full polypeptide chain with name, 3-letter, and 1-letter codes |
| Protein | Composition stats, charge, hydrophobicity analysis |
| Database lookup | Queries UniProt BLAST for real protein matches |

---

## Sample Sequences to Test

**DNA (non-template strand):**
```
ATGAAACCCGGGTTTTAA
```

**DNA (template strand):**
```
TTAATGCGATCGATCGAT
```

**RNA:**
```
AUGAAACCCGGGUUUUAA
```

**Invalid:**
```
ATGXYZ123
```

---

## Hosting on the Internet (Phase E equivalent)

### Render (Free)
1. Push folder to GitHub
2. Go to [render.com](https://render.com) → New Web Service
3. Connect repo, set:
   - **Build:** `pip install -r requirements.txt`
   - **Start:** `gunicorn app:app`
4. Done — get a public URL

---

## Marking Scheme Coverage

| Component | Marks | Implemented |
|---|---|---|
| Sequence Input | 10 | Text area, file upload, drag & drop |
| Sequence Detection | 15 | DNA/RNA/invalid detection + plain-English explanation |
| Transcription | 15 | Both strand types, mRNA shown, explanation |
| Translation | 20 | Full codon table, start/stop codons, per-codon display |
| Amino Acids | 15 | Full polypeptide chain with names, 3-letter, 1-letter codes |
| Protein & Database | 20 | Composition stats + UniProt BLAST lookup |
| Overall Quality | 5 | Clean UI, all explanations present |
| **Total** | **100** | |
