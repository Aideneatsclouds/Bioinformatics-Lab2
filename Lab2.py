import streamlit as st
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Configure Entrez Email
Entrez.email = "ythinaharan@gmail.com"

# --- Page Configuration ---
st.set_page_config(page_title="Bioinformatics Lab Suite", layout="wide")


# ==========================================
# LAB 1 FUNCTIONS: Retrieval & Analysis
# ==========================================
def retrieve_data(protein_id):
    """Retrieves a protein sequence record from NCBI Entrez using UniProt ID."""
    try:
        handle = Entrez.efetch(
            db="protein", id=protein_id, rettype="fasta", retmode="text"
        )
        record = SeqIO.read(handle, "fasta")
        return record
    except Exception:
        return None


def get_basic_analysis(sequence):
    """Analyzes sequence length, amino acid composition, molecular weight, and pI."""
    seq_analysis = ProteinAnalysis(str(sequence))

    length = len(sequence)
    aa_composition = seq_analysis.count_amino_acids()
    molecular_weight = seq_analysis.molecular_weight()
    isoelectric_point = seq_analysis.isoelectric_point()

    return length, aa_composition, molecular_weight, isoelectric_point


# ==========================================
# LAB 2 FUNCTIONS: Pairwise Alignment (From Scratch)
# ==========================================
def global_alignment(seq1, seq2, match_score, mismatch_score, gap_penalty):
    """Needleman-Wunsch Global Alignment Algorithm implemented from scratch."""
    m, n = len(seq1), len(seq2)

    # Initialize scoring matrix and traceback pointer matrix
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    # Initialize first column and first row
    for i in range(1, m + 1):
        score_matrix[i][0] = i * gap_penalty
        traceback[i][0] = 'U'  # Up
    for j in range(1, n + 1):
        score_matrix[0][j] = j * gap_penalty
        traceback[0][j] = 'L'  # Left

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            max_score = max(match, delete, insert)
            score_matrix[i][j] = max_score

            # Record direction for backtrace
            if max_score == match:
                traceback[i][j] = 'D'  # Diagonal
            elif max_score == delete:
                traceback[i][j] = 'U'  # Up
            else:
                traceback[i][j] = 'L'  # Left

    # Backtrace from bottom-right cell (m, n)
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    backtrace_path = [(i, j)]

    while i > 0 or j > 0:
        direction = traceback[i][j]
        if direction == 'D':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif direction == 'L':
            aligned_seq1.append('-')
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

    # Initialize scoring matrix (all zeros) and traceback matrix
    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    traceback = [[None for _ in range(n + 1)] for _ in range(m + 1)]

    max_score = 0
    max_pos = (0, 0)

    # Fill the scoring matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix[i - 1][j - 1] + (
                match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
            )
            delete = score_matrix[i - 1][j] + gap_penalty
            insert = score_matrix[i][j - 1] + gap_penalty

            # Local alignment resets negative scores to 0
            cell_score = max(0, match, delete, insert)
            score_matrix[i][j] = cell_score

            if cell_score == 0:
                traceback[i][j] = None
            elif cell_score == match:
                traceback[i][j] = 'D'
            elif cell_score == delete:
                traceback[i][j] = 'U'
            else:
                traceback[i][j] = 'L'

            # Track highest score in the entire matrix
            if cell_score >= max_score:
                max_score = cell_score
                max_pos = (i, j)

    # Backtrace starting from the highest-scoring cell until a 0 is reached
    aligned_seq1, aligned_seq2 = [], []
    i, j = max_pos
    backtrace_path = [(i, j)]

    while i > 0 and j > 0 and score_matrix[i][j] > 0:
        direction = traceback[i][j]
        if direction == 'D':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif direction == 'U':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif direction == 'L':
            aligned_seq1.append('-')
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


def style_alignment_matrix(matrix, seq1, seq2, path):
    """Creates a styled Pandas DataFrame highlighting the backtrace path."""
    # Add space padding for the 0th row and column
    row_labels = ["-"] + list(seq1)
    col_labels = ["-"] + list(seq2)

    # Create unique multi-index or formatted labels to prevent duplicate index styling issues
    df = pd.DataFrame(matrix)
    path_set = set(path)

    def highlight_cells(data):
        styles = pd.DataFrame("", index=data.index, columns=data.columns)
        for r, c in path_set:
            styles.iloc[r, c] = (
                "background-color: #2e7d32; color: white; font-weight: bold;"
            )
        return styles

    # Apply styling using integer coordinates before renaming headers
    styled_df = df.style.apply(highlight_cells, axis=None)

    # Set readable row/column headers (with index number to guarantee uniqueness)
    styled_df.index = [f"{idx}: {char}" for idx, char in enumerate(row_labels)]
    styled_df.columns = [f"{idx}: {char}" for idx, char in enumerate(col_labels)]

    return styled_df


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

# Initialize session state to optionally pass sequences from Lab 1 to Lab 2
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
                st.session_state.retrieved_seq = str(record.seq)

                # Display basic record information
                st.subheader("1. Protein Record Details")
                st.write(f"**Name / ID:** `{record.name}`")
                st.write(f"**Description:** {record.description}")

                with st.expander("View Full FASTA Sequence"):
                    st.code(str(record.seq), language="text")

                # Perform sequence analysis
                length, aa_comp, mol_weight, iso_point = get_basic_analysis(
                    record.seq
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

    # Option to load sequence from Lab 1 if available
    if st.session_state.retrieved_seq:
        st.info(
            f"A sequence from Lab 1 ({len(st.session_state.retrieved_seq)} aa) is"
            " available in memory."
        )
        if st.button("Use First 15 Residues of Lab 1 Sequence as Sequence 1"):
            st.session_state.seq1_input = st.session_state.retrieved_seq[:15]

    # --- Alignment Configuration ---
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

    # --- Scoring Metrics ---
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
            # Execute chosen algorithm
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

            # Display Alignment Score
            st.metric(label="Optimal Alignment Score", value=score)

            # Display Aligned Sequences
            st.markdown("**Aligned Sequences:**")
            alignment_visual = format_alignment_display(aln1, aln2)
            st.code(alignment_visual, language="text")

            # Display Backtrace Path Coordinates
            with st.expander("View Backtrace Coordinate Path (Row, Col)"):
                st.write(" -> ".join([f"({r}, {c})" for r, c in path]))

            # Display Alignment Matrix with Highlighted Backtrace
            st.markdown("**Alignment Matrix & Highlighted Backtrace Path:**")
            st.caption(
                "Green cells indicate the optimal backtrace path. Row/Column"
                " headers show `Index: Character`."
            )

            styled_matrix = style_alignment_matrix(matrix, seq1, seq2, path)
            st.dataframe(styled_matrix, use_container_width=True)