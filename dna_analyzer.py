"""
ATGo — DNA Sequence Analyzer
=============================
From Start Codon to Discovery.

A bioinformatics tool that analyzes DNA sequences.

Features:
1. GC Content - measures DNA stability
2. Reverse Complement - generates the opposite DNA strand
3. Open Reading Frames (ORFs) - finds potential protein-coding regions
4. Mutation Detection - compares two sequences for substitutions,
   insertions, deletions, and frameshifts
"""


# ============================================================
# PART 1: FASTA PARSER
# ============================================================
# FASTA is a standard text format for DNA sequences.
# It looks like this:
#   >Sequence_Name
#   ATGCGATCGATCG...
#
# The first line starts with ">" and is a label (header).
# Everything after is the DNA sequence.
# Our parser reads this format and extracts just the DNA letters.

def parse_fasta(fasta_text):
    """
    Takes a string in FASTA format and returns:
    - header: the name/description line (without the >)
    - sequence: the DNA letters joined into one string (uppercase)

    Example:
        Input:  ">My_Gene\nATGCGATCG\nATCGATCG"
        Output: ("My_Gene", "ATGCGATCGATCGATCG")
    """
    lines = fasta_text.strip().split("\n")  # Split text into lines
    header = ""
    sequence_parts = []  # We'll collect sequence lines here

    for line in lines:
        line = line.strip()  # Remove extra spaces/tabs
        if line.startswith(">"):
            # This is the header line — store it without the ">"
            header = line[1:].strip()
        elif line:
            # This is a sequence line — add it to our collection
            sequence_parts.append(line.upper())

    # Join all sequence parts into one continuous string
    sequence = "".join(sequence_parts)

    # Validate: make sure it only contains A, T, G, C
    valid_bases = set("ATGC")
    for base in sequence:
        if base not in valid_bases:
            print(f"  Warning: Found unexpected character '{base}' in sequence.")
            print(f"  Only A, T, G, C are valid DNA bases.")
            # Remove invalid characters
            sequence = "".join(b for b in sequence if b in valid_bases)
            print(f"  Cleaned sequence length: {len(sequence)} bases")
            break

    return header, sequence


# ============================================================
# PART 2: GC CONTENT ANALYZER
# ============================================================
# GC content = percentage of bases that are G or C.
#
# WHY IT MATTERS (real biology):
# - G-C pairs form 3 hydrogen bonds (stronger)
# - A-T pairs form 2 hydrogen bonds (weaker)
# - More GC = more stable DNA = harder to separate strands
# - Used in PCR primer design, species classification, etc.

def calculate_gc_content(sequence):
    """
    Calculates the percentage of G and C bases in a DNA sequence.

    Formula: GC% = (count of G + count of C) / total bases * 100

    Example:
        Input:  "ATGCGC"
        G count = 2, C count = 2, total = 6
        GC% = (2 + 2) / 6 * 100 = 66.67%
    """
    if len(sequence) == 0:
        return 0.0

    g_count = sequence.count("G")  # Count how many G's
    c_count = sequence.count("C")  # Count how many C's
    total = len(sequence)           # Total number of bases

    gc_percentage = ((g_count + c_count) / total) * 100
    return round(gc_percentage, 2)  # Round to 2 decimal places


def interpret_gc(gc_percentage):
    """
    Gives a biological interpretation of the GC percentage.
    These thresholds are based on real scientific observations:

    - Thermophilic organisms (heat-loving) have very high GC
      because their DNA needs to stay stable at high temperatures.
    - The human genome is about 41% GC.
    - Some parasitic bacteria have very low GC (~25%).
    """
    if gc_percentage > 65:
        return ("Very high GC — extremely stable DNA. "
                "Typical of thermophilic organisms (heat-loving bacteria) "
                "that live in hot springs or volcanic vents.")
    elif gc_percentage > 55:
        return ("High GC — stable DNA. "
                "Common in soil bacteria like Streptomyces "
                "(these bacteria produce many antibiotics!).")
    elif gc_percentage > 45:
        return ("Moderate GC — normal range for many organisms. "
                "Many bacteria and some eukaryotes fall here.")
    elif gc_percentage > 35:
        return ("Moderate-low GC — normal for many animals. "
                "The human genome is about 41% GC.")
    else:
        return ("Low GC — AT-rich DNA, less thermally stable. "
                "Seen in some parasitic organisms like Mycoplasma "
                "and Plasmodium (malaria parasite).")


def display_gc_content(sequence):
    """Runs GC analysis and prints a nice report."""
    print("\n" + "=" * 55)
    print("  GC CONTENT ANALYSIS")
    print("=" * 55)

    total = len(sequence)
    g_count = sequence.count("G")
    c_count = sequence.count("C")
    a_count = sequence.count("A")
    t_count = sequence.count("T")
    gc = calculate_gc_content(sequence)

    print(f"\n  Total bases:    {total}")
    print(f"  A (Adenine):    {a_count}  ({round(a_count/total*100, 1)}%)")
    print(f"  T (Thymine):    {t_count}  ({round(t_count/total*100, 1)}%)")
    print(f"  G (Guanine):    {g_count}  ({round(g_count/total*100, 1)}%)")
    print(f"  C (Cytosine):   {c_count}  ({round(c_count/total*100, 1)}%)")
    print(f"\n  GC Content:     {gc}%")
    print(f"  AT Content:     {round(100 - gc, 2)}%")
    print(f"\n  Interpretation: {interpret_gc(gc)}")


# ============================================================
# PART 3: REVERSE COMPLEMENT
# ============================================================
# DNA is double-stranded and the two strands run in opposite
# directions (antiparallel). The complement rules are:
#   A <-> T    (they pair together)
#   C <-> G    (they pair together)
#
# To get the reverse complement:
#   Step 1: Replace each base with its complement
#   Step 2: Reverse the entire string
#
# WHY IT MATTERS:
# - Genes can be on either strand of DNA
# - To analyze the other strand, you need the reverse complement
# - Essential for primer design in PCR

def reverse_complement(sequence):
    """
    Returns the reverse complement of a DNA sequence.

    Example:
        Input:    5'-ATGCGA-3'
        Step 1 (complement): TACGCT
        Step 2 (reverse):    TCGCAT
        Output:   3'-TCGCAT-5'  (which is read 5'->3' as TCGCAT)
    """
    # Dictionary mapping each base to its complement
    complement_map = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G"
    }

    # Step 1: Replace each base with its complement
    complemented = ""
    for base in sequence:
        complemented += complement_map.get(base, base)

    # Step 2: Reverse the string
    # In Python, [::-1] reverses a string
    rev_comp = complemented[::-1]

    return rev_comp


def display_reverse_complement(sequence):
    """Shows the reverse complement with a visual diagram."""
    print("\n" + "=" * 55)
    print("  REVERSE COMPLEMENT")
    print("=" * 55)

    rev_comp = reverse_complement(sequence)

    # Show a short preview (first 50 bases) for readability
    display_len = min(len(sequence), 50)
    original_display = sequence[:display_len]
    revcomp_display = rev_comp[:display_len]
    suffix = "..." if len(sequence) > 50 else ""

    print(f"\n  Original (5'→3'):     {original_display}{suffix}")
    print(f"  Rev. Complement (5'→3'): {revcomp_display}{suffix}")
    print(f"\n  Full length: {len(rev_comp)} bases")
    print(f"\n  Full reverse complement:")
    # Print in lines of 60 characters (standard format)
    for i in range(0, len(rev_comp), 60):
        print(f"  {rev_comp[i:i+60]}")


# ============================================================
# PART 4: OPEN READING FRAME (ORF) FINDER
# ============================================================
# An ORF is a stretch of DNA that COULD code for a protein.
#
# Rules:
#   - Starts with ATG (the start codon = Methionine)
#   - Ends with TAA, TAG, or TGA (stop codons)
#   - Is read in groups of 3 bases (codons)
#
# DNA can be read in 3 different "reading frames":
#   Frame 1: ATG|CGA|TCG|...  (start at position 0)
#   Frame 2:  A|TGC|GAT|CG... (start at position 1)
#   Frame 3:  AT|GCG|ATC|G... (start at position 2)
#
# WHY IT MATTERS:
# - ORFs help predict where genes are in a genome
# - If an ORF is long, it's likely a real gene
# - Used in genome annotation (labeling what each part of DNA does)

# Codon table: translates 3-letter DNA codons to amino acids
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}


def translate_to_protein(dna_sequence):
    """
    Translates a DNA sequence into a protein (amino acid) sequence.
    Reads in groups of 3 (codons) and looks up each one in the codon table.
    Stops when it hits a stop codon (* symbol).
    """
    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, "?")
        if amino_acid == "*":
            break  # Stop codon reached
        protein += amino_acid
    return protein


def find_orfs(sequence, min_length=30):
    """
    Finds all Open Reading Frames in a DNA sequence.

    Searches all 3 reading frames for patterns: ATG ... (stop codon)

    Parameters:
        sequence: DNA string (A, T, G, C)
        min_length: minimum ORF length in bases (default 30 = 10 amino acids)
                    Short ORFs are usually not real genes.

    Returns a list of dictionaries, each containing:
        - frame: which reading frame (1, 2, or 3)
        - start: start position in the sequence
        - end: end position
        - length: length in bases
        - protein: translated amino acid sequence
    """
    orfs = []

    # Check all 3 reading frames
    for frame in range(3):
        i = frame  # Start position for this frame

        while i < len(sequence) - 2:
            codon = sequence[i:i+3]

            if codon == "ATG":  # Found a start codon!
                # Now scan forward looking for a stop codon
                orf_start = i
                j = i + 3  # Move to next codon

                while j < len(sequence) - 2:
                    next_codon = sequence[j:j+3]
                    if next_codon in STOP_CODONS:
                        # Found a stop codon — this is the end of the ORF
                        orf_end = j + 3  # Include the stop codon
                        orf_sequence = sequence[orf_start:orf_end]

                        if len(orf_sequence) >= min_length:
                            protein = translate_to_protein(orf_sequence)
                            orfs.append({
                                "frame": frame + 1,
                                "start": orf_start + 1,  # 1-indexed for biology
                                "end": orf_end,
                                "length": len(orf_sequence),
                                "codons": len(orf_sequence) // 3,
                                "protein": protein,
                                "dna": orf_sequence
                            })
                        i = j + 3  # Continue searching after this ORF
                        break
                    j += 3
                else:
                    # No stop codon found — move on
                    i += 3
                    continue
            else:
                i += 3  # Not a start codon, move to next codon

    # Sort by length (longest first — longer ORFs are more likely real genes)
    orfs.sort(key=lambda x: x["length"], reverse=True)
    return orfs


def display_orfs(sequence):
    """Finds and displays all ORFs with protein translations."""
    print("\n" + "=" * 55)
    print("  OPEN READING FRAME (ORF) FINDER")
    print("=" * 55)

    # Also check the reverse complement strand
    rev_comp = reverse_complement(sequence)

    print("\n  Searching forward strand (+)...")
    forward_orfs = find_orfs(sequence)

    print(f"  Searching reverse strand (-)...")
    reverse_orfs = find_orfs(rev_comp)

    # Label them
    for orf in forward_orfs:
        orf["strand"] = "+"
    for orf in reverse_orfs:
        orf["strand"] = "-"

    all_orfs = forward_orfs + reverse_orfs
    all_orfs.sort(key=lambda x: x["length"], reverse=True)

    if not all_orfs:
        print("\n  No ORFs found (minimum length: 30 bases / 10 amino acids).")
        print("  Try a longer sequence or lower the minimum length.")
        return

    print(f"\n  Found {len(all_orfs)} ORF(s):\n")

    for i, orf in enumerate(all_orfs, 1):
        print(f"  --- ORF #{i} ---")
        print(f"  Strand:        {orf['strand']} ({'forward' if orf['strand'] == '+' else 'reverse complement'})")
        print(f"  Reading Frame: {orf['frame']}")
        print(f"  Position:      {orf['start']} to {orf['end']}")
        print(f"  Length:         {orf['length']} bases ({orf['codons']} codons)")
        print(f"  Protein ({len(orf['protein'])} aa): {orf['protein']}")

        # Show first 60 bases of the DNA
        dna_preview = orf["dna"][:60]
        suffix = "..." if len(orf["dna"]) > 60 else ""
        print(f"  DNA:           {dna_preview}{suffix}")
        print()


# ============================================================
# PART 5: MUTATION DETECTION
# ============================================================
# Compares a REFERENCE sequence (normal/healthy) to a
# TEST sequence (patient/sample) and finds:
#
#   1. Substitutions: one base replaced by another (e.g., A → G)
#      - Most common type of mutation
#      - Can be silent, missense, or nonsense
#
#   2. Insertions: extra bases added into the sequence
#      - If not a multiple of 3, causes a FRAMESHIFT
#
#   3. Deletions: bases removed from the sequence
#      - If not a multiple of 3, causes a FRAMESHIFT
#
#   4. Frameshifts: insertions or deletions that shift the
#      reading frame, usually destroying the protein
#
# WHY IT MATTERS:
# - Mutations cause diseases (cancer, sickle cell, cystic fibrosis)
# - Detecting mutations helps with diagnosis and treatment
# - Used in personalized medicine and genetic counseling

def detect_mutations(reference, test):
    """
    Compares two DNA sequences and finds all mutations.

    Uses position-by-position comparison for equal-length sequences,
    and detects length differences as insertions/deletions.

    Parameters:
        reference: the "normal" DNA sequence
        test: the "sample" DNA sequence to compare

    Returns a dictionary with:
        - substitutions: list of single-base changes
        - insertions: extra bases in the test sequence
        - deletions: missing bases in the test sequence
        - frameshifts: insertions/deletions not divisible by 3
        - summary: overall statistics
    """
    mutations = {
        "substitutions": [],
        "insertions": [],
        "deletions": [],
        "frameshifts": [],
    }

    ref_len = len(reference)
    test_len = len(test)

    # --- Detect insertions and deletions based on length ---
    if test_len > ref_len:
        # Test is longer → insertion detected
        extra_bases = test_len - ref_len
        insertion = {
            "type": "insertion",
            "position": ref_len + 1,
            "bases_added": extra_bases,
            "inserted_sequence": test[ref_len:],
            "is_frameshift": (extra_bases % 3) != 0
        }
        mutations["insertions"].append(insertion)

        if insertion["is_frameshift"]:
            mutations["frameshifts"].append({
                "type": "frameshift (insertion)",
                "bases": extra_bases,
                "detail": f"{extra_bases} base(s) inserted — NOT a multiple of 3, "
                          f"shifts reading frame!"
            })

    elif test_len < ref_len:
        # Test is shorter → deletion detected
        missing_bases = ref_len - test_len
        deletion = {
            "type": "deletion",
            "position": test_len + 1,
            "bases_deleted": missing_bases,
            "deleted_sequence": reference[test_len:],
            "is_frameshift": (missing_bases % 3) != 0
        }
        mutations["deletions"].append(deletion)

        if deletion["is_frameshift"]:
            mutations["frameshifts"].append({
                "type": "frameshift (deletion)",
                "bases": missing_bases,
                "detail": f"{missing_bases} base(s) deleted — NOT a multiple of 3, "
                          f"shifts reading frame!"
            })

    # --- Detect substitutions (compare matching positions) ---
    compare_length = min(ref_len, test_len)
    for i in range(compare_length):
        if reference[i] != test[i]:
            mutations["substitutions"].append({
                "position": i + 1,  # 1-indexed for biology
                "reference_base": reference[i],
                "test_base": test[i],
                "change": f"{reference[i]} → {test[i]}"
            })

    # --- Build summary ---
    total_mutations = (len(mutations["substitutions"]) +
                       len(mutations["insertions"]) +
                       len(mutations["deletions"]))

    mutations["summary"] = {
        "total_mutations": total_mutations,
        "substitutions": len(mutations["substitutions"]),
        "insertions": len(mutations["insertions"]),
        "deletions": len(mutations["deletions"]),
        "frameshifts": len(mutations["frameshifts"]),
        "ref_length": ref_len,
        "test_length": test_len,
        "mutation_rate": round((len(mutations["substitutions"]) / compare_length * 100), 2)
        if compare_length > 0 else 0
    }

    return mutations


def display_mutations(reference, test):
    """Runs mutation detection and prints a detailed report."""
    print("\n" + "=" * 55)
    print("  MUTATION DETECTION")
    print("=" * 55)

    results = detect_mutations(reference, test)
    summary = results["summary"]

    # Show sequence comparison
    display_len = min(len(reference), len(test), 50)
    print(f"\n  Reference: {reference[:display_len]}{'...' if len(reference) > 50 else ''}")
    print(f"  Test:      {test[:display_len]}{'...' if len(test) > 50 else ''}")

    # Show alignment markers (| for match, X for mismatch)
    markers = ""
    for i in range(display_len):
        if reference[i] == test[i]:
            markers += "|"
        else:
            markers += "X"
    print(f"  Match:     {markers}{'...' if min(len(reference), len(test)) > 50 else ''}")

    # Summary
    print(f"\n  --- Summary ---")
    print(f"  Reference length: {summary['ref_length']} bases")
    print(f"  Test length:      {summary['test_length']} bases")
    print(f"  Substitutions:    {summary['substitutions']}")
    print(f"  Insertions:       {summary['insertions']}")
    print(f"  Deletions:        {summary['deletions']}")
    print(f"  Frameshifts:      {summary['frameshifts']}")
    print(f"  Mutation rate:    {summary['mutation_rate']}%")

    # Detailed substitutions
    if results["substitutions"]:
        print(f"\n  --- Substitutions ---")
        for sub in results["substitutions"]:
            print(f"  Position {sub['position']}: {sub['change']}")

    # Detailed insertions
    if results["insertions"]:
        print(f"\n  --- Insertions ---")
        for ins in results["insertions"]:
            preview = ins["inserted_sequence"][:30]
            suffix = "..." if len(ins["inserted_sequence"]) > 30 else ""
            print(f"  At position {ins['position']}: "
                  f"+{ins['bases_added']} base(s) inserted: {preview}{suffix}")
            if ins["is_frameshift"]:
                print(f"  ⚠ FRAMESHIFT: This insertion shifts the reading frame!")

    # Detailed deletions
    if results["deletions"]:
        print(f"\n  --- Deletions ---")
        for dele in results["deletions"]:
            preview = dele["deleted_sequence"][:30]
            suffix = "..." if len(dele["deleted_sequence"]) > 30 else ""
            print(f"  At position {dele['position']}: "
                  f"-{dele['bases_deleted']} base(s) deleted: {preview}{suffix}")
            if dele["is_frameshift"]:
                print(f"  ⚠ FRAMESHIFT: This deletion shifts the reading frame!")

    # Frameshift warning
    if results["frameshifts"]:
        print(f"\n  ⚠ FRAMESHIFT WARNING ⚠")
        print(f"  Frameshifts change the reading frame of ALL downstream codons.")
        print(f"  This usually produces a completely different (nonfunctional) protein.")
        for fs in results["frameshifts"]:
            print(f"  → {fs['detail']}")

    if summary["total_mutations"] == 0:
        print(f"\n  ✓ No mutations detected — sequences are identical!")


# ============================================================
# PART 6: MAIN MENU
# ============================================================
# This ties everything together into an interactive program.
# The user picks what they want to do from a menu.

def get_sequence_input():
    """
    Asks the user to input a DNA sequence.
    Accepts either:
      - Raw DNA (e.g., ATGCGATCG)
      - FASTA format (e.g., >Gene_Name\nATGCGATCG)
    """
    print("\n  Enter your DNA sequence (FASTA format or raw sequence).")
    print("  For multi-line input, type 'END' on a new line when done:")

    lines = []
    while True:
        line = input("  ")
        if line.strip().upper() == "END":
            break
        lines.append(line)

    full_input = "\n".join(lines)

    # Check if it's FASTA format (starts with >)
    if full_input.strip().startswith(">"):
        header, sequence = parse_fasta(full_input)
        print(f"\n  Parsed FASTA sequence: '{header}'")
    else:
        header = "User Input"
        sequence = full_input.replace("\n", "").replace(" ", "").upper()
        # Clean invalid characters
        valid = set("ATGC")
        sequence = "".join(b for b in sequence if b in valid)

    print(f"  Sequence length: {len(sequence)} bases")
    return header, sequence


def main():
    """Main program loop — shows menu and runs selected analysis."""
    print("\n" + "=" * 55)
    print("  🧬 ATGo — DNA Sequence Analyzer 🧬")
    print("  From Start Codon to Discovery")
    print("=" * 55)
    print("  A tool for analyzing DNA sequences using Python.")
    print("  Data source: NCBI GenBank (ncbi.nlm.nih.gov/genbank)")

    while True:
        print("\n" + "-" * 55)
        print("  MAIN MENU")
        print("-" * 55)
        print("  1. Analyze GC Content")
        print("  2. Get Reverse Complement")
        print("  3. Find Open Reading Frames (ORFs)")
        print("  4. Detect Mutations (compare 2 sequences)")
        print("  5. Run ALL analyses on one sequence")
        print("  6. Run demo with sample data")
        print("  7. Exit")
        print("-" * 55)

        choice = input("  Enter your choice (1-7): ").strip()

        if choice == "1":
            _, seq = get_sequence_input()
            if seq:
                display_gc_content(seq)

        elif choice == "2":
            _, seq = get_sequence_input()
            if seq:
                display_reverse_complement(seq)

        elif choice == "3":
            _, seq = get_sequence_input()
            if seq:
                display_orfs(seq)

        elif choice == "4":
            print("\n  --- Reference Sequence (normal/healthy) ---")
            _, ref_seq = get_sequence_input()
            print("\n  --- Test Sequence (patient/sample) ---")
            _, test_seq = get_sequence_input()
            if ref_seq and test_seq:
                display_mutations(ref_seq, test_seq)

        elif choice == "5":
            _, seq = get_sequence_input()
            if seq:
                display_gc_content(seq)
                display_reverse_complement(seq)
                display_orfs(seq)

        elif choice == "6":
            run_demo()

        elif choice == "7":
            print("\n  Thanks for using DNA Sequence Analyzer!")
            print("  Keep exploring bioinformatics! 🧬\n")
            break

        else:
            print("  Invalid choice. Please enter 1-7.")


def run_demo():
    """
    Runs a demonstration with a real-world sample sequence.
    This is a short segment inspired by the human BRCA1 gene,
    which is involved in breast cancer when mutated.
    """
    print("\n" + "=" * 55)
    print("  DEMO MODE — Sample Analysis")
    print("=" * 55)
    print("  Using a sample sequence for demonstration.")

    # Sample sequence (a short made-up but realistic sequence)
    sample_fasta = """>Sample_Gene_Demo
ATGGCTGACATCGATCGTAGCTAGCATGCGAATCGATCGATATCGATCGTAGCTA
GCATGCGATCGATCGATCGATCGATCGATCGATCGTAGCTAGCATGCGATCGATC
GATCGATCGATCGATCGTAGCTAGCATGCGATCGATCGATCGATCGTAAATGCCC
GGGAAATTTCCCGGGATCGATCGATCGATCGATCGATCGTAGCTAGCTGA"""

    header, sequence = parse_fasta(sample_fasta)
    print(f"  Sequence: {header}")
    print(f"  Length: {len(sequence)} bases\n")

    # Run all analyses
    display_gc_content(sequence)
    display_reverse_complement(sequence)
    display_orfs(sequence)

    # Demo mutation detection
    print("\n  --- Mutation Detection Demo ---")
    reference = "ATGGCTGACATCGATCGTAGCTAGC"
    #                    ^         ^  ^
    test_mutated = "ATGGCTGACATCAATCGTGGCTAGCAA"  # substitutions + insertion
    print(f"  Reference: {reference}")
    print(f"  Mutated:   {test_mutated}")
    display_mutations(reference, test_mutated)


# ============================================================
# RUN THE PROGRAM
# ============================================================
# This is the entry point. When you run:
#   python dna_analyzer.py
# Python starts here and calls main()

if __name__ == "__main__":
    main()
