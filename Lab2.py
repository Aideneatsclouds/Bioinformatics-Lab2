import re
import streamlit as st
from Bio import Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Configure Entrez Email
Entrez.email = "ythinaharan@gmail.com"

# --- Page Configuration ---
st.set_page_config(page_title="Bioinformatics Lab Suite", layout="wide")


# ==========================================
# LAB 1 FUNCTIONS: Retrieval & Analysis
# ==========================================
def parse_fasta_text(fasta_text):
    """Parses raw FASTA text in pure Python without importing Bio.SeqIO."""
    lines = [line.strip() for line in fasta_text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        return None

    header = lines[0][1:].strip()
    record_id = header.split()[0] if header else "Unknown"
    sequence = "".join(line.replace(" ", "") for line in lines[1:] if not line.startswith(">")).upper()

    if not sequence:
        return None

    return {
        "name": record_id,
        "description": header,
        "seq": sequence,
    }


def retrieve_data(protein_id):
    """Retrieves a protein sequence record from NCBI Entrez using UniProt ID."""
    try:
        handle = Entrez.efetch(
            db="protein", id=protein_id, rettype="fasta", retmode="text"
        )
        fasta_text = handle.read()
        handle.close()
        if isinstance(fasta_text, bytes):
            fasta_text = fasta_text.decode("utf-8", errors="replace")
        return parse_fasta_text(fasta_text)
    except Exception:
        return None


# Standard amino acid molecular weights (average Da, water-adjusted) and pKa values
AA_WEIGHTS = {
    "A": 89.09, "R": 174.20, "N": 132.12, "D": 133.10, "C": 121.15,
    "E": 147.13, "Q": 146.15, "G": 75.07, "H": 155.16, "I": 131.17,
    "L": 131.17, "K": 146.19, "M": 149.21, "F": 165.19, "P": 115.13,
    "S": 105.09, "T": 119.12, "W": 204.23, "Y": 181.19, "V": 117.15,
}
H2O_WEIGHT = 18.015


def estimate_pi_pure_python(seq_str):
    """Estimates isoelectric point (pI) using standard pKa values if ProtParam is unavailable."""
    counts = {aa: seq_str.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}
    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2.0
        # Positive charges: N-term (~9.69), K (~10.5), R (~12.48), H (~6.0)
        pos = (
            1.0 / (1.0 + 10.0 ** (mid - 9.69))
            + counts["K"] / (1.0 + 10.0 ** (mid - 10.5))
            + counts["R"] / (1.0 + 10.0 ** (mid - 12.48))
            + counts["H"] / (1.0 + 10.0 ** (mid - 6.0))
        )
        # Negative charges: C-term (~2.34), D (~3.65), E (~4.25), C (~8.18), Y (~10.07)
        neg = (
            1.0 / (1.0 + 10.0 ** (2.34 - mid))
            + counts["D"] / (1.0 + 10.0 ** (3.65 - mid))
            + counts["E"] / (1.0 + 10.0 ** (4.25 - mid))
            + counts["C"] / (1.0 + 10.0 ** (8.18 - mid))
            + counts["Y"] / (1.0 + 10.0 ** (10.07 - mid))
        )
        if pos - neg > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2.0


def get_basic_analysis(sequence):
    """Analyzes sequence length, amino acid composition, molecular weight, and pI."""
    seq_str = str(sequence).upper()
    length = len(seq_str)

    try:
        seq_analysis = ProteinAnalysis(seq_str)
        aa_composition = seq_analysis.count_amino_acids()
        molecular_weight = seq_analysis.molecular_weight()
        isoelectric_point = seq_analysis.isoelectric_point()
    except Exception:
        cleaned = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", seq_str)
        aa_composition = {aa: cleaned.count(aa) for aa in "ACDEFGHIKLMNPQRSTVWY"}
        total_aa = sum(aa_composition.values())
        molecular_weight = (
            sum(count * AA_WEIGHTS[aa] for aa, count in aa_composition.items())
            - max(0, total_aa - 1) * H2O_WEIGHT
        )
        isoelectric_point = estimate_pi_pure_python(cleaned)

    return length, aa_composition, molecular_weight, isoelectric_point


# ==========================================
# LAB 2 FUNCTIONS: Pairwise Alignment (From Scratch)
# ==========================================
def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
    """Needleman-Wunsch Global Alignment Algorithm implemented from scratch."""
    m, n = len(seq1), len(seq2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback[i][0] = "U"
    for j in range(1, n + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback[0][j] = "L"

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            max_score = max(match, delete, insert)
            score_matrix[i][j] = max_score

            if max_score == match:
                traceback[i][j] = "D"
            elif max_score == delete:
                traceback[i][j] = "U"
            else:
                traceback[i][j] = "L"

    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    backtrace_path = [(i, j)]

    while i > 0 or j > 0:
        direction = traceback[i][j]
        if direction == "D":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif direction == "L":
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        backtrace_path.append((i, j))

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))
    final_score = score_matrix[m][n]

    return final_score, aligned_seq1, aligned_seq2, score_matrix, backtrace_path


def local_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
    """Smith-Waterman Local Alignment Algorithm implemented from scratch."""
    m, n = len(seq1), len(seq2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            cell_score = max(0, match, delete, insert)
            score_matrix[i][j] = cell_score

            if cell_score == 0:
                traceback[i][j] = None
            elif cell_score == match:
                traceback[i][j] = "D"
            elif cell_score == delete:
                traceback[i][j] = "U"
            else:
                traceback[i][j] = "L"

            if cell_score >= max_score:
                max_score = cell_score
                max_pos = (i, j)

    aligned_seq1, aligned_seq2 = [], []
    i, j = max_pos
    backtrace_path = [(i, j)]

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        direction = traceback[i][j]
        if direction == "D":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == "U":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1
        elif direction == "L":
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            break
        backtrace_path.append((i, j))

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return max_score, aligned_seq1, aligned_seq2, score_matrix, backtrace_path


def format_alignment_display(aligned1, aligned2):
    """Generates a visual match/mismatch midline between two aligned sequences."""
    midline = []
    for a, b in zip(aligned1, aligned2):
        if a == b:
            midline.append("|")
        elif a == "-" or b == "-":
            midline.append(" ")
        else:
            midline.append(".")
    return f"Seq 1: {aligned1}\n       {''.join(midline)}\nSeq 2: {aligned2}"


def render_alignment_matrix_html(matrix, seq1, seq2, path):
    """Renders the scoring matrix as an HTML table with backtrace cells highlighted."""
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)
    path_set = set(path)

    html = [
        "<div style='overflow-x: auto;'>",
        "<table style='border-collapse: collapse; text-align: center; font-family: monospace;'>",
    ]

    html.append("<tr>")
    html.append(
        "<th style='border: 1px solid #555; padding: 6px 10px; background: #222; color: #fff;'></th>"
    )
    for j, col_char in enumerate(col_labels):
        html.append(
            f"<th style='border: 1px solid #555; padding: 6px 10px; background: #222; color: #fff;'>{j}: {col_char}</th>"
        )
    html.append("</tr>")

    for i, row_char in enumerate(row_labels):
        html.append("<tr>")
        html.append(
            f"<th style='border: 1px solid #555; padding: 6px 10px; background: #222; color: #fff;'>{i}: {row_char}</th>"
        )
        for j, val in enumerate(matrix[i]):
            if (i, j) in path_set:
                cell_style = (
                    "border: 1px solid #555; padding: 6px 10px; "
                    "background-color: #2e7d32; color: #ffffff; font-weight: bold;"
                )
            else:
                cell_style = "border: 1px solid #555; padding: 6px 10px;"
            html.append(f"<td style='{cell_style}'>{val}</td>")
        html.append("</tr>")

    html.append("</table></div>")
    return "".join(html)


# ==========================================
# SIDEBAR NAVIGATION & HEADER
# ==========================================
st.sidebar.title("Navigation")
st.sidebar.markdown("**Course:** Bioinformatics Lab")
st.sidebar.markdown("**Student:** Yaaswiny Abbiramavally")
st.sidebar.divider()

lab_selection = st.sidebar.radio(
    "Select Lab Module:",
    ["Lab 1: Sequence Retrieval & Analysis", "Lab 2: Pairwise Sequence Alignment"],
)

if "retrieved_seq" not in st.session_state:
    st.session_state.retrieved_seq = ""


# ==========================================
# VIEW 1: LAB 1
# ==========================================
if lab_selection == "Lab 1: Sequence Retrieval & Analysis":
    st.title("Lab 1 - Protein Sequence Retrieval & Analysis")
    st.markdown(
        "**Course:** Bioinformatics Lab | **Student:** Yaaswiny Abbiramavally"
    )

    protein_id = st.text_input("Enter UniProt ID", value="P04637")
    retrieve_btn = st.button("Retrieve and Analyze")

    if retrieve_btn:
        if protein_id.strip() != "":
            with st.spinner("Fetching sequence data from database..."):
                record = retrieve_data(protein_id.strip())

            if record:
                st.success("Successfully retrieved protein record!")
                st.session_state.retrieved_seq = record["seq"]

                st.subheader("1. Protein Record Details")
                st.write(f"**Name / ID:** `{record['name']}`")
                st.write(f"**Description:** {record['description']}")

                with st.expander("View Full FASTA Sequence"):
                    st.code(record["seq"], language="text")

                length, aa_comp, mol_weight, iso_point = get_basic_analysis(
                    record["seq"]
                )

                st.subheader("2. Basic Sequence Analysis")
                col1, col2 = st.columns(2)
                with col1:
                    st.metric(label="Sequence Length", value=f"{length} aa")
                    st.metric(
                        label="Isoelectric Point (pI)", value=f"{iso_point:.2f}"
                    )

                with col2:
                    st.metric(
                        label="Molecular Weight", value=f"{mol_weight:,.2f} Da"
                    )

                st.markdown("### Amino Acid Composition")
                st.json(aa_comp)
            else:
                st.error(
                    "Failed to retrieve record. Please check the UniProt ID and"
                    " try again."
                )
        else:
            st.warning("Please enter a valid UniProt ID.")


# ==========================================
# VIEW 2: LAB 2
# ==========================================
elif lab_selection == "Lab 2: Pairwise Sequence Alignment":
    st.title("Lab 2 - Pairwise Sequence Alignment")
    st.markdown(
        "Custom implementation of **Needleman-Wunsch (Global)** and"
        " **Smith-Waterman (Local)** algorithms."
    )

    if st.session_state.retrieved_seq:
        st.info(
            f"A sequence from Lab 1 ({len(st.session_state.retrieved_seq)} aa) is"
            " available in memory."
        )
        if st.button("Use First 15 Residues of Lab 1 Sequence as Sequence 1"):
            st.session_state.seq1_input = st.session_state.retrieved_seq[:15]

    st.subheader("1. Alignment Settings & Input Sequences")

    alignment_type = st.radio(
        "Select Alignment Type:",
        ["Global Alignment (Needleman-Wunsch)", "Local Alignment (Smith-Waterman)"],
        horizontal=True,
    )

    col_seq1, col_seq2 = st.columns(2)
    with col_seq1:
        seq1 = st.text_input(
            "Sequence 1 (Vertical / Rows)",
            value=st.session_state.get("seq1_input", "GCATGCG"),
        ).strip().upper()
    with col_seq2:
        seq2 = st.text_input(
            "Sequence 2 (Horizontal / Columns)",
            value="GATTACA",
        ).strip().upper()

    st.markdown("**Scoring Parameters:**")
    col_m, col_mm, col_g = st.columns(3)
    with col_m:
        match_score = st.number_input("Match Score", value=1, step=1)
    with col_mm:
        mismatch_score = st.number_input("Mismatch Score", value=-1, step=1)
    with col_g:
        gap_penalty = st.number_input("Gap Penalty", value=-1, step=1)

    align_btn = st.button("Run Sequence Alignment", type="primary")

    if align_btn:
        if not seq1 or not seq2:
            st.error("Please enter both Sequence 1 and Sequence 2.")
        else:
            if "Global" in alignment_type:
                score, aln1, aln2, matrix, path = global_alignment(
                    seq1, seq2, match_score, mismatch_score, gap_penalty
                )
            else:
                score, aln1, aln2, matrix, path = local_alignment(
                    seq1, seq2, match_score, mismatch_score, gap_penalty
                )

            st.divider()
            st.subheader("2. Alignment Results")

            st.metric(label="Optimal Alignment Score", value=score)

            st.markdown("**Aligned Sequences:**")
            alignment_visual = format_alignment_display(aln1, aln2)
            st.code(alignment_visual, language="text")

            with st.expander("View Backtrace Coordinate Path (Row, Col)"):
                st.write(" -> ".join([f"({r}, {c})" for r, c in path]))

            st.markdown("**Alignment Matrix & Highlighted Backtrace Path:**")
            st.caption(
                "Green cells indicate the optimal backtrace path. Row/Column"
                " headers show `Index: Character`."
            )

            matrix_html = render_alignment_matrix_html(matrix, seq1, seq2, path)
            st.markdown(matrix_html, unsafe_allow_html=True)