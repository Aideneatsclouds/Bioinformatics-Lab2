import re
import streamlit as st
from Bio import Entrez
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Configure Entrez Email
Entrez.email = "ythinaharan@gmail.com"

# --- Page Configuration ---
st.set_page_config(page_title="Bioinformatics Lab Suite", layout="wide")


# ==========================================
# CUSTOM CSS: Pastel Pink & Matcha Green Theme
# ==========================================
def apply_custom_theme():
    st.markdown(
        """
        <style>
        /* Global Font & Pastel Gradient Background */
        html, body, [class*="css"], .stApp, p, label, span, div, input, button, textarea {
            font-family: 'Times New Roman', Times, serif !important;
        }

        .stApp {
            background: linear-gradient(135deg, #fce4ec 0%, #f8bbd0 35%, #e8f5e9 70%, #c8e6c9 100%) !important;
            color: #3d2c35 !important;
        }

        /* Sidebar Styling */
        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #f48fb1 0%, #a5d6a7 100%) !important;
            border-right: 2px solid rgba(255, 255, 255, 0.7);
            border-top-right-radius: 24px;
            border-bottom-right-radius: 24px;
        }

        [data-testid="stSidebar"] * {
            color: #ffffff !important;
            text-shadow: 0 0 8px rgba(255, 255, 255, 0.85), 0 1px 3px rgba(90, 45, 65, 0.45);
        }

        /* Glowing White Headers */
        h1, h2, h3, [data-testid="stHeading"] * {
            font-family: 'Times New Roman', Times, serif !important;
            color: #ffffff !important;
            text-shadow:
                0 0 8px rgba(255, 255, 255, 0.95),
                0 0 18px rgba(244, 143, 177, 0.9),
                0 2px 6px rgba(136, 74, 102, 0.55) !important;
            letter-spacing: 0.4px;
        }

        h1 {
            background: linear-gradient(90deg, rgba(244, 143, 177, 0.75), rgba(165, 214, 167, 0.75));
            padding: 14px 22px !important;
            border-radius: 20px !important;
            border: 1.5px solid rgba(255, 255, 255, 0.85);
            box-shadow: 0 6px 20px rgba(244, 143, 177, 0.25);
        }

        h2, h3 {
            background: rgba(240, 168, 194, 0.6);
            padding: 8px 16px !important;
            border-radius: 16px !important;
            border: 1px solid rgba(255, 255, 255, 0.75);
            display: inline-block;
            margin-top: 8px !important;
        }

        /* Rounded Cards for Metrics, Alerts, and Expanders */
        [data-testid="stMetric"] {
            background: rgba(255, 255, 255, 0.78) !important;
            border: 2px solid #f48fb1 !important;
            border-radius: 20px !important;
            padding: 16px 20px !important;
            box-shadow: 0 6px 16px rgba(165, 214, 167, 0.35);
        }

        [data-testid="stMetricLabel"] * {
            color: #6d4c5d !important;
            font-size: 1.05rem !important;
        }

        [data-testid="stMetricValue"] * {
            color: #4e7d55 !important;
            font-weight: bold !important;
        }

        /* Inputs & Number Boxes */
        div[data-baseweb="input"] > div,
        div[data-baseweb="base-input"] {
            background-color: rgba(255, 255, 255, 0.88) !important;
            border: 1.5px solid #a5d6a7 !important;
            border-radius: 16px !important;
            color: #3d2c35 !important;
        }

        input {
            color: #3d2c35 !important;
            font-size: 1.05rem !important;
        }

        /* Rounded Pastel Buttons */
        .stButton > button {
            background: linear-gradient(90deg, #f48fb1 0%, #a5d6a7 100%) !important;
            color: #ffffff !important;
            font-weight: bold !important;
            font-size: 1.05rem !important;
            border: 1.5px solid #ffffff !important;
            border-radius: 22px !important;
            padding: 10px 24px !important;
            text-shadow: 0 0 8px rgba(255, 255, 255, 0.9), 0 1px 2px rgba(80, 40, 60, 0.4);
            box-shadow: 0 4px 14px rgba(244, 143, 177, 0.4);
            transition: all 0.25s ease-in-out;
        }

        .stButton > button:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 20px rgba(129, 199, 132, 0.55);
            border-color: #fce4ec !important;
        }

        /* Code Blocks & Expanders */
        pre, code {
            font-family: 'Times New Roman', Times, serif !important;
            background-color: rgba(255, 255, 255, 0.85) !important;
            color: #3d2c35 !important;
            border-radius: 16px !important;
            border: 1.5px solid #f8bbd0 !important;
            font-size: 1.05rem !important;
        }

        [data-testid="stExpander"] {
            background-color: rgba(255, 255, 255, 0.72) !important;
            border: 1.5px solid #c8e6c9 !important;
            border-radius: 18px !important;
        }

        [data-testid="stAlert"] {
            border-radius: 18px !important;
            border: 1.5px solid rgba(255, 255, 255, 0.85) !important;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


apply_custom_theme()


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
    sequence = "".join(
        line.replace(" ", "") for line in lines[1:] if not line.startswith(">")
    ).upper()

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
        pos = (
            1.0 / (1.0 + 10.0 ** (mid - 9.69))
            + counts["K"] / (1.0 + 10.0 ** (mid - 10.5))
            + counts["R"] / (1.0 + 10.0 ** (mid - 12.48))
            + counts["H"] / (1.0 + 10.0 ** (mid - 6.0))
        )
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
    """Renders the scoring matrix in a pastel pink & matcha theme with glowing backtrace cells."""
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)
    path_set = set(path)

    html = [
        "<div style='overflow-x: auto; padding: 12px; background: rgba(255,255,255,0.75); "
        "border-radius: 20px; border: 2px solid #f48fb1; box-shadow: 0 6px 18px rgba(244,143,177,0.25);'>",
        "<table style='border-collapse: separate; border-spacing: 4px; text-align: center; "
        "font-family: \"Times New Roman\", Times, serif; margin: auto;'>",
    ]

    # Header row
    html.append("<tr>")
    html.append(
        "<th style='padding: 8px 12px; border-radius: 12px; background: #f48fb1; color: #ffffff; "
        "text-shadow: 0 0 8px rgba(255,255,255,0.9); border: 1px solid #ffffff;'></th>"
    )
    for j, col_char in enumerate(col_labels):
        html.append(
            f"<th style='padding: 8px 12px; border-radius: 12px; background: #f48fb1; color: #ffffff; "
            f"text-shadow: 0 0 8px rgba(255,255,255,0.9); border: 1px solid #ffffff;'>{j}: {col_char}</th>"
        )
    html.append("</tr>")

    # Data rows
    for i, row_char in enumerate(row_labels):
        html.append("<tr>")
        html.append(
            f"<th style='padding: 8px 12px; border-radius: 12px; background: #f48fb1; color: #ffffff; "
            f"text-shadow: 0 0 8px rgba(255,255,255,0.9); border: 1px solid #ffffff;'>{i}: {row_char}</th>"
        )
        for j, val in enumerate(matrix[i]):
            if (i, j) in path_set:
                cell_style = (
                    "padding: 8px 12px; border-radius: 12px; "
                    "background: linear-gradient(135deg, #81c784, #66bb6a); "
                    "color: #ffffff; font-weight: bold; "
                    "text-shadow: 0 0 8px rgba(255,255,255,0.95); "
                    "border: 1.5px solid #ffffff; "
                    "box-shadow: 0 0 10px rgba(102, 187, 106, 0.7);"
                )
            else:
                cell_style = (
                    "padding: 8px 12px; border-radius: 12px; "
                    "background-color: #fff5f8; color: #4a3540; "
                    "border: 1px solid #f8bbd0;"
                )
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

                st.subheader("3. Amino Acid Composition")
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

    align_btn = st.button("Run Sequence Alignment")

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

            st.subheader("3. Alignment Matrix & Backtrace Path")
            st.caption(
                "Glowing matcha-green cells indicate the optimal backtrace path."
                " Row/Column headers show Index: Character."
            )

            matrix_html = render_alignment_matrix_html(matrix, seq1, seq2, path)
            st.markdown(matrix_html, unsafe_allow_html=True)