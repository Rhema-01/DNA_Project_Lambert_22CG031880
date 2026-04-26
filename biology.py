"""
biology.py — Core DNA/RNA analysis engine
CSC 442 - Project 2
"""

# ── Codon table (mRNA codons → amino acid) ────────────────────────────────────

CODON_TABLE = {
    # Phenylalanine
    "UUU": ("Phenylalanine", "Phe", "F"),
    "UUC": ("Phenylalanine", "Phe", "F"),
    # Leucine
    "UUA": ("Leucine",       "Leu", "L"),
    "UUG": ("Leucine",       "Leu", "L"),
    "CUU": ("Leucine",       "Leu", "L"),
    "CUC": ("Leucine",       "Leu", "L"),
    "CUA": ("Leucine",       "Leu", "L"),
    "CUG": ("Leucine",       "Leu", "L"),
    # Isoleucine
    "AUU": ("Isoleucine",    "Ile", "I"),
    "AUC": ("Isoleucine",    "Ile", "I"),
    "AUA": ("Isoleucine",    "Ile", "I"),
    # Methionine / START
    "AUG": ("Methionine",    "Met", "M"),
    # Valine
    "GUU": ("Valine",        "Val", "V"),
    "GUC": ("Valine",        "Val", "V"),
    "GUA": ("Valine",        "Val", "V"),
    "GUG": ("Valine",        "Val", "V"),
    # Serine
    "UCU": ("Serine",        "Ser", "S"),
    "UCC": ("Serine",        "Ser", "S"),
    "UCA": ("Serine",        "Ser", "S"),
    "UCG": ("Serine",        "Ser", "S"),
    "AGU": ("Serine",        "Ser", "S"),
    "AGC": ("Serine",        "Ser", "S"),
    # Proline
    "CCU": ("Proline",       "Pro", "P"),
    "CCC": ("Proline",       "Pro", "P"),
    "CCA": ("Proline",       "Pro", "P"),
    "CCG": ("Proline",       "Pro", "P"),
    # Threonine
    "ACU": ("Threonine",     "Thr", "T"),
    "ACC": ("Threonine",     "Thr", "T"),
    "ACA": ("Threonine",     "Thr", "T"),
    "ACG": ("Threonine",     "Thr", "T"),
    # Alanine
    "GCU": ("Alanine",       "Ala", "A"),
    "GCC": ("Alanine",       "Ala", "A"),
    "GCA": ("Alanine",       "Ala", "A"),
    "GCG": ("Alanine",       "Ala", "A"),
    # Tyrosine
    "UAU": ("Tyrosine",      "Tyr", "Y"),
    "UAC": ("Tyrosine",      "Tyr", "Y"),
    # Stop codons
    "UAA": ("STOP",          "***", "*"),
    "UAG": ("STOP",          "***", "*"),
    "UGA": ("STOP",          "***", "*"),
    # Histidine
    "CAU": ("Histidine",     "His", "H"),
    "CAC": ("Histidine",     "His", "H"),
    # Glutamine
    "CAA": ("Glutamine",     "Gln", "Q"),
    "CAG": ("Glutamine",     "Gln", "Q"),
    # Asparagine
    "AAU": ("Asparagine",    "Asn", "N"),
    "AAC": ("Asparagine",    "Asn", "N"),
    # Lysine
    "AAA": ("Lysine",        "Lys", "K"),
    "AAG": ("Lysine",        "Lys", "K"),
    # Aspartate
    "GAU": ("Aspartate",     "Asp", "D"),
    "GAC": ("Aspartate",     "Asp", "D"),
    # Glutamate
    "GAA": ("Glutamate",     "Glu", "E"),
    "GAG": ("Glutamate",     "Glu", "E"),
    # Cysteine
    "UGU": ("Cysteine",      "Cys", "C"),
    "UGC": ("Cysteine",      "Cys", "C"),
    # Tryptophan
    "UGG": ("Tryptophan",    "Trp", "W"),
    # Arginine
    "CGU": ("Arginine",      "Arg", "R"),
    "CGC": ("Arginine",      "Arg", "R"),
    "CGA": ("Arginine",      "Arg", "R"),
    "CGG": ("Arginine",      "Arg", "R"),
    "AGA": ("Arginine",      "Arg", "R"),
    "AGG": ("Arginine",      "Arg", "R"),
    # Glycine
    "GGU": ("Glycine",       "Gly", "G"),
    "GGC": ("Glycine",       "Gly", "G"),
    "GGA": ("Glycine",       "Gly", "G"),
    "GGG": ("Glycine",       "Gly", "G"),
}

DNA_BASES  = set("ATCG")
RNA_BASES  = set("AUCG")
STOP_CODONS = {"UAA", "UAG", "UGA"}
START_CODON = "AUG"

# ── Sequence detection ────────────────────────────────────────────────────────

def clean_sequence(raw: str) -> str:
    """Remove whitespace, FASTA headers, newlines; uppercase."""
    lines = raw.strip().splitlines()
    cleaned = []
    for line in lines:
        if line.startswith(">"):   # FASTA header
            continue
        cleaned.append(line.strip().upper().replace(" ", ""))
    return "".join(cleaned)


def detect_sequence(seq: str) -> dict:
    """
    Returns dict with keys: type ('DNA'|'RNA'|'INVALID'), invalid_chars, explanation.
    """
    chars = set(seq)
    has_T = "T" in chars
    has_U = "U" in chars
    invalid = chars - (DNA_BASES | {"U"})   # everything outside ATCGU

    if invalid:
        return {
            "type": "INVALID",
            "invalid_chars": sorted(invalid),
            "explanation": (
                f"The sequence contains characters that do not belong to either DNA or RNA: "
                f"{', '.join(sorted(invalid))}. "
                "Valid DNA uses only A, T, C, G. Valid RNA uses only A, U, C, G."
            )
        }

    if has_T and has_U:
        return {
            "type": "INVALID",
            "invalid_chars": ["T", "U"],
            "explanation": (
                "The sequence contains both T (Thymine) and U (Uracil). "
                "DNA uses T but never U; RNA uses U but never T. "
                "A real sequence cannot contain both."
            )
        }

    if has_U:
        seq_type = "RNA"
        explanation = (
            "The sequence was identified as RNA because it contains the base U (Uracil). "
            "Uracil is found only in RNA — DNA uses Thymine (T) instead. "
            "All bases in the sequence (A, U, C, G) are valid RNA bases."
        )
    else:
        seq_type = "DNA"
        explanation = (
            "The sequence was identified as DNA because it contains only the bases A, T, C, and G — "
            "the four bases that make up DNA. "
            "The presence of T (Thymine) in particular confirms this is DNA, "
            "since RNA uses U (Uracil) in place of T."
        )

    return {"type": seq_type, "invalid_chars": [], "explanation": explanation}


# ── Transcription ─────────────────────────────────────────────────────────────

def transcribe(seq: str, seq_type: str, strand_type: str = "non-template") -> dict:
    """
    Produce mRNA from a DNA or RNA sequence.
    strand_type: 'non-template' (coding strand) or 'template'
    """
    if seq_type == "RNA":
        mrna = seq  # already mRNA
        explanation = (
            "The input is already RNA (specifically mRNA), so no transcription step is needed. "
            "The sequence is used directly for translation."
        )
        return {"mrna": mrna, "input_seq": seq, "explanation": explanation}

    if strand_type == "non-template":
        # Non-template (coding) strand: same sequence as mRNA, just replace T→U
        mrna = seq.replace("T", "U")
        explanation = (
            f"You provided the non-template strand (also called the coding or sense strand). "
            f"Transcription works by replacing each Thymine (T) with Uracil (U) to produce mRNA. "
            f"No complementary step is needed for this strand type. "
            f"Input: {seq[:40]}{'...' if len(seq)>40 else ''} → "
            f"mRNA:  {mrna[:40]}{'...' if len(mrna)>40 else ''}"
        )
    else:
        # Template strand: take complement, then replace T→U
        complement_map = {"A": "U", "T": "A", "C": "G", "G": "C"}
        mrna = "".join(complement_map.get(b, b) for b in seq)
        explanation = (
            f"You provided the template strand (also called the antisense strand). "
            f"During transcription, RNA polymerase reads the template strand and builds the "
            f"complementary mRNA strand: A pairs with U, T pairs with A, C pairs with G, G pairs with C. "
            f"Input (template): {seq[:40]}{'...' if len(seq)>40 else ''} → "
            f"mRNA: {mrna[:40]}{'...' if len(mrna)>40 else ''}"
        )

    return {"mrna": mrna, "input_seq": seq, "explanation": explanation}


# ── Translation ───────────────────────────────────────────────────────────────

def translate(mrna: str) -> dict:
    """
    Translate mRNA into codons and amino acids.
    Starts at first AUG, stops at stop codon or end of sequence.
    """
    # Find start codon
    start_idx = mrna.find(START_CODON)
    if start_idx == -1:
        return {
            "codons": [],
            "amino_acids": [],
            "protein_seq": "",
            "found_start": False,
            "explanation": (
                "No start codon (AUG) was found in the mRNA sequence. "
                "Translation cannot begin without a start codon. "
                "In real biology, no protein would be produced from this sequence."
            )
        }

    coding_region = mrna[start_idx:]
    codons = []
    amino_acids = []
    protein_letters = []

    for i in range(0, len(coding_region) - 2, 3):
        codon = coding_region[i:i+3]
        if len(codon) < 3:
            break

        info = CODON_TABLE.get(codon, ("Unknown", "???", "?"))
        name, three_letter, one_letter = info
        codons.append({
            "codon": codon,
            "name": name,
            "three": three_letter,
            "one": one_letter,
            "is_start": (codon == START_CODON and i == 0),
            "is_stop":  (codon in STOP_CODONS)
        })

        if codon in STOP_CODONS:
            amino_acids.append({"name": "STOP", "three": "***", "one": "*", "codon": codon})
            break
        else:
            amino_acids.append({"name": name, "three": three_letter, "one": one_letter, "codon": codon})
            protein_letters.append(one_letter)

    protein_seq = "".join(protein_letters)

    explanation = (
        f"Translation begins at the first AUG codon (start codon) found at position {start_idx + 1} "
        f"of the mRNA. The ribosome reads the mRNA three bases at a time — each group of three is called "
        f"a codon. Each codon tells the ribosome which amino acid to add to the growing chain. "
        f"Translation stops when a stop codon (UAA, UAG, or UGA) is reached. "
        f"A total of {len(amino_acids)} codons were read, producing a chain of "
        f"{len(protein_letters)} amino acid(s)."
    )

    return {
        "codons": codons,
        "amino_acids": amino_acids,
        "protein_seq": protein_seq,
        "found_start": True,
        "start_idx": start_idx,
        "explanation": explanation
    }


# ── Protein characterisation ──────────────────────────────────────────────────

def characterise_protein(protein_seq: str, amino_acids: list) -> dict:
    """Basic composition and property analysis of the protein."""
    if not protein_seq:
        return {"composition": {}, "properties": {}, "explanation": "No protein sequence to characterise."}

    # Amino acid composition
    composition = {}
    for aa in amino_acids:
        if aa["name"] == "STOP":
            continue
        name = aa["name"]
        composition[name] = composition.get(name, 0) + 1

    # Simple property estimation
    # Charged residues
    positive = sum(1 for aa in amino_acids if aa["one"] in "KRH")
    negative = sum(1 for aa in amino_acids if aa["one"] in "DE")
    hydrophobic = sum(1 for aa in amino_acids if aa["one"] in "AVILMFYW")
    polar = sum(1 for aa in amino_acids if aa["one"] in "STNQ")
    total = len(protein_seq)

    charge = positive - negative
    if charge > 0:
        charge_desc = f"Net positive (basic) — {positive} positively charged, {negative} negatively charged residues"
    elif charge < 0:
        charge_desc = f"Net negative (acidic) — {positive} positively charged, {negative} negatively charged residues"
    else:
        charge_desc = f"Neutral — {positive} positively charged, {negative} negatively charged residues"

    explanation = (
        f"A protein is the final product of the gene expression pipeline. It is formed when the "
        f"polypeptide chain folds into a specific 3D shape determined by its amino acid sequence. "
        f"The protein produced here is {total} amino acid(s) long. "
        f"It has a {charge_desc.lower()}. "
        f"{round(hydrophobic/total*100) if total else 0}% of its residues are hydrophobic "
        f"(water-repelling), and {round(polar/total*100) if total else 0}% are polar (water-attracting). "
        f"These properties influence how the protein folds and what function it might perform."
    )

    return {
        "composition": composition,
        "length": total,
        "positive_charged": positive,
        "negative_charged": negative,
        "hydrophobic": hydrophobic,
        "polar": polar,
        "net_charge": charge,
        "charge_desc": charge_desc,
        "explanation": explanation,
    }
