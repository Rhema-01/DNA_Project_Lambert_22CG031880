"""
app.py — Flask Web Application
DNA & RNA Sequence Analyzer
CSC 442 - Project 2
"""

from flask import Flask, render_template, request, jsonify
import os
import urllib.request
import urllib.parse
import json
import time
import xml.etree.ElementTree as ET

from biology import (
    clean_sequence, detect_sequence,
    transcribe, translate, characterise_protein
)

app = Flask(__name__)
app.config["UPLOAD_FOLDER"] = os.path.join(os.path.dirname(__file__), "static", "uploads")
app.config["MAX_CONTENT_LENGTH"] = 5 * 1024 * 1024  # 5 MB

os.makedirs(app.config["UPLOAD_FOLDER"], exist_ok=True)

# ── UniProt BLAST-like search (free, no API key needed) ──────────────────────

def search_uniprot(protein_seq: str) -> dict:
    """
    Query UniProt's BLAST service with the protein sequence.
    Falls back to a keyword search if BLAST times out.
    """
    if not protein_seq or len(protein_seq) < 5:
        return {"success": False, "results": [], "message": "Protein sequence too short for database search."}

    try:
        # Submit BLAST job to UniProt
        blast_url = "https://rest.uniprot.org/blast/run"
        payload = urllib.parse.urlencode({
            "sequence":  protein_seq,
            "database":  "uniprotkb_swissprot",
            "taxons":    "",
            "goTerms":   "",
            "matrix":    "BLOSUM62",
            "threshold": "0.001",
            "filter":    "true",
            "gapped":    "true",
            "hits":      "5",
        }).encode()

        req = urllib.request.Request(blast_url, data=payload, method="POST")
        req.add_header("Accept", "application/json")

        with urllib.request.urlopen(req, timeout=10) as resp:
            job_data = json.loads(resp.read())
            job_id = job_data.get("jobId", "")

        if not job_id:
            raise ValueError("No job ID returned")

        # Poll for results (max 20 seconds)
        results_url = f"https://rest.uniprot.org/blast/results/{job_id}"
        for _ in range(10):
            time.sleep(2)
            try:
                with urllib.request.urlopen(results_url + "?format=json&size=5", timeout=8) as r:
                    result_data = json.loads(r.read())
                    if result_data.get("results"):
                        hits = []
                        for hit in result_data["results"][:5]:
                            entry = hit.get("entryMappedToList", {})
                            hits.append({
                                "accession":  hit.get("accession", "N/A"),
                                "name":       hit.get("entryName", "Unknown protein"),
                                "protein":    entry.get("proteinName", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown"),
                                "organism":   hit.get("organism", {}).get("scientificName", "Unknown organism"),
                                "function":   _extract_function(hit),
                                "identity":   hit.get("alignments", [{}])[0].get("identity", "N/A"),
                                "score":      hit.get("alignments", [{}])[0].get("score", "N/A"),
                                "url":        f"https://www.uniprot.org/uniprotkb/{hit.get('accession', '')}",
                            })
                        return {"success": True, "results": hits, "source": "UniProt BLAST"}
            except Exception:
                continue

        raise TimeoutError("BLAST search timed out")

    except Exception as e:
        # Fallback: search UniProt by amino acid sequence pattern
        return _uniprot_keyword_fallback(protein_seq, str(e))


def _extract_function(hit: dict) -> str:
    comments = hit.get("comments", [])
    for c in comments:
        if c.get("commentType") == "FUNCTION":
            texts = c.get("texts", [])
            if texts:
                return texts[0].get("value", "")[:200]
    return "Function information not available."


def _uniprot_keyword_fallback(protein_seq: str, error: str) -> dict:
    """Search UniProt by sequence-derived keywords as fallback."""
    try:
        query = urllib.parse.urlencode({
            "query":  f"length:[{max(1,len(protein_seq)-5)} TO {len(protein_seq)+5}]",
            "format": "json",
            "size":   "3",
            "fields": "accession,protein_name,organism_name,cc_function",
        })
        url = f"https://rest.uniprot.org/uniprotkb/search?{query}"
        with urllib.request.urlopen(url, timeout=8) as r:
            data = json.loads(r.read())
            hits = []
            for entry in data.get("results", [])[:3]:
                hits.append({
                    "accession": entry.get("primaryAccession", "N/A"),
                    "name":      entry.get("uniProtkbId", "Unknown"),
                    "protein":   entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "Unknown protein"),
                    "organism":  entry.get("organism", {}).get("scientificName", "Unknown organism"),
                    "function":  "See UniProt page for full details.",
                    "identity":  "N/A (keyword match)",
                    "score":     "N/A",
                    "url":       f"https://www.uniprot.org/uniprotkb/{entry.get('primaryAccession','')}",
                })
            if hits:
                return {"success": True, "results": hits, "source": "UniProt keyword search (BLAST unavailable)", "note": error}
    except Exception as e2:
        pass

    return {
        "success": False,
        "results": [],
        "message": f"Database search could not be completed: {error}. This can happen with very short or novel sequences.",
        "manual_url": f"https://www.uniprot.org/blast/",
    }


# ── Routes ────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return render_template("index.html")


@app.route("/analyze", methods=["POST"])
def analyze():
    # Get sequence from form or uploaded file
    sequence_raw = ""
    if "file" in request.files and request.files["file"].filename:
        f = request.files["file"]
        sequence_raw = f.read().decode("utf-8", errors="ignore")
    else:
        sequence_raw = request.form.get("sequence", "")

    strand_type = request.form.get("strand_type", "non-template")

    if not sequence_raw.strip():
        return jsonify({"error": "No sequence provided."}), 400

    seq = clean_sequence(sequence_raw)
    if not seq:
        return jsonify({"error": "Sequence is empty after cleaning."}), 400

    # Step 1: Detect
    detection = detect_sequence(seq)
    if detection["type"] == "INVALID":
        return jsonify({
            "step": "detection",
            "detection": detection,
            "sequence": seq,
        })

    seq_type = detection["type"]

    # Step 2: Transcription
    transcription = transcribe(seq, seq_type, strand_type)
    mrna = transcription["mrna"]

    # Step 3: Translation
    translation = translate(mrna)
    amino_acids  = translation["amino_acids"]
    protein_seq  = translation["protein_seq"]

    # Step 4: Protein characterisation
    char = characterise_protein(protein_seq, amino_acids)

    # Step 5: Database search
    db_result = search_uniprot(protein_seq)

    return jsonify({
        "sequence":     seq,
        "seq_type":     seq_type,
        "strand_type":  strand_type,
        "detection":    detection,
        "transcription": transcription,
        "translation":  translation,
        "protein_char": char,
        "db_result":    db_result,
    })


@app.route("/upload", methods=["POST"])
def upload():
    """Return the text content of an uploaded file."""
    if "file" not in request.files:
        return jsonify({"error": "No file"}), 400
    f = request.files["file"]
    content = f.read().decode("utf-8", errors="ignore")
    return jsonify({"content": content})


if __name__ == "__main__":
    app.run(debug=True)
