# ============================================
# 1. IMPORT LIBRARIES
# ============================================

import streamlit as st


# ============================================
# 2. PAGE CONFIG
# ============================================

st.set_page_config(
    page_title="DNA Cloning Calculator | RNMolBioTools",
    page_icon="🧬",
    layout="wide"
)


# ============================================
# 3. HEADER
# ============================================

st.markdown(
"""
<div style="
background: linear-gradient(90deg,#2C3E50,#4CA1AF);
padding:16px;
border-radius:8px;
display:flex;
align-items:center;
gap:12px;
color:white;
box-shadow:0px 2px 6px rgba(0,0,0,0.2);
">

<div style="font-size:34px;">🧬</div>

<div>
<div style="font-size:20px;font-weight:600;">
DNA Cloning Calculator
</div>

<div style="font-size:13px;opacity:0.9;">
RNMolBioTools
</div>
</div>

</div>
""",
unsafe_allow_html=True
)

st.write("")


# ============================================
# 4. NAVIGATION
# ============================================

if "page" not in st.session_state:
    st.session_state.page = "cloning"

col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    if st.button("🧬 DNA Cloning", use_container_width=True):
        st.session_state.page = "cloning"

with col2:
    if st.button("🔢 Copy Number", use_container_width=True):
        st.session_state.page = "copy"

with col3:
    if st.button("⚖ Mass Converter", use_container_width=True):
        st.session_state.page = "conversion"

with col4:
    if st.button("🧾 Sequence Mass", use_container_width=True):
        st.session_state.page = "sequence"

with col5:
    if st.button("🧪 Sequence Tools", use_container_width=True):
        st.session_state.page = "tools"

st.markdown("---")


# ============================================
# 5. DNA CLONING CALCULATOR
# ============================================

if st.session_state.page == "cloning":

    st.markdown("### DNA Cloning Calculator")

    vector_kb = st.number_input("Vector DNA size (kb)", min_value=0.0)
    insert_kb = st.number_input("Insert DNA size (kb)", min_value=0.0)
    ratio = st.text_input("Vector : Insert ratio; adjust accordingly", "1:4")
    total_ng = st.number_input("Total DNA (ng); adjust accordingly", value=100.0)

    if st.button("Calculate Cloning Mix"):

        r_vec, r_ins = ratio.split(":")
        r_vec, r_ins = float(r_vec), float(r_ins)

        vector_bp = vector_kb * 1000
        insert_bp = insert_kb * 1000

        ratio_val = r_ins / r_vec

        vector_ng = total_ng / (1 + ratio_val * (insert_bp / vector_bp))
        insert_ng = total_ng - vector_ng

        vector_pmol = (vector_ng * 1000) / (vector_bp * 650)
        insert_pmol = (insert_ng * 1000) / (insert_bp * 650)

        st.subheader("Results")

        st.write(f"Vector: {vector_ng:.2f} ng ({vector_pmol:.4f} pmol)")
        st.write(f"Insert: {insert_ng:.2f} ng ({insert_pmol:.4f} pmol)")

    st.info("""
Typical ligation ratios

Sticky end ligations: **1 : 3**  
Blunt end ligations: **1 : 4 – 6**  
Large fragment ligations: **1 : 1**
""")

    st.caption("ⓘ Calculation assumptions")

    with st.expander("How the ligation calculation works"):

        st.markdown("""
DNA ligation efficiency depends on **molar ratios** rather than mass.

Average molecular weight of dsDNA:

**650–660 g/mol per bp**

Conversion used:

pmol = (ng × 1000) / (bp × 650)
""")


# ============================================
# 6. DNA COPY NUMBER CALCULATOR
# ============================================

elif st.session_state.page == "copy":

    st.markdown("### dsDNA Copy Number Calculator")

    dna_kb = st.number_input("DNA length (kb)", min_value=0.0)
    dna_ng = st.number_input("DNA amount (ng)", min_value=0.0)

    if st.button("Calculate Copy Number"):

        dna_bp = dna_kb * 1000
        molar_mass = dna_bp * 650

        moles = (dna_ng * 1e-9) / molar_mass
        pmol = moles * 1e12
        copies = moles * 6.022e23

        st.subheader("Results")

        st.write(f"Molar Mass: {molar_mass:.2e} g/mol")
        st.write(f"Moles: {moles:.2e} mol")
        st.write(f"pmol: {pmol:.4f}")
        st.write(f"DNA Copies: {copies:.2e}")

    st.caption("ⓘ Calculation assumptions")

    with st.expander("How DNA copy number is calculated"):

        st.markdown("""
Average molecular weight:

**650 g/mol per base pair**

Copies are calculated using **Avogadro's number**:

6.022 × 10²³ molecules/mol
""")


# ============================================
# 7. DNA / RNA MASS CONVERSION
# ============================================

elif st.session_state.page == "conversion":

    st.markdown("### DNA/RNA Mass Unit Converter")

    na_type = st.selectbox("Molecule type", ["dsDNA", "ssDNA", "RNA"])

    for key in ["g","mg","ug","ng","pg"]:
        if key not in st.session_state:
            st.session_state[key] = 0.0

    def convert_from_g():
        g = st.session_state.g
        st.session_state.mg = g * 1e3
        st.session_state.ug = g * 1e6
        st.session_state.ng = g * 1e9
        st.session_state.pg = g * 1e12

    def convert_from_mg():
        mg = st.session_state.mg
        st.session_state.g = mg / 1e3
        st.session_state.ug = mg * 1e3
        st.session_state.ng = mg * 1e6
        st.session_state.pg = mg * 1e9

    def convert_from_ug():
        ug = st.session_state.ug
        st.session_state.g = ug / 1e6
        st.session_state.mg = ug / 1e3
        st.session_state.ng = ug * 1e3
        st.session_state.pg = ug * 1e6

    def convert_from_ng():
        ng = st.session_state.ng
        st.session_state.g = ng / 1e9
        st.session_state.mg = ng / 1e6
        st.session_state.ug = ng / 1e3
        st.session_state.pg = ng * 1e3

    def convert_from_pg():
        pg = st.session_state.pg
        st.session_state.g = pg / 1e12
        st.session_state.mg = pg / 1e9
        st.session_state.ug = pg / 1e6
        st.session_state.ng = pg / 1e3

    col1, col2, col3, col4, col5 = st.columns(5)

    col1.number_input("g", key="g", on_change=convert_from_g)
    col2.number_input("mg", key="mg", on_change=convert_from_mg)
    col3.number_input("µg", key="ug", on_change=convert_from_ug)
    col4.number_input("ng", key="ng", on_change=convert_from_ng)
    col5.number_input("pg", key="pg", on_change=convert_from_pg)

    st.write("Enter a value in any field to convert automatically.")


# ============================================
# 8. DNA SEQUENCE MASS CALCULATOR
# ============================================

elif st.session_state.page == "sequence":

    st.markdown("### DNA Sequence Mass Calculator")

    seq_type = st.selectbox("Sequence Type", ["dsDNA", "ssDNA"])

    sequence = st.text_area("Paste DNA sequence (A, T, G, C only)", height=150)

    mass_ng = st.number_input("Optional: DNA amount (ng)", min_value=0.0, value=0.0)

    if sequence:

        seq = sequence.upper().replace("\n","").replace(" ","")

        valid_bases = set("ATGC")

        if not set(seq).issubset(valid_bases):

            st.error("Invalid input: sequence must contain only A, T, G, C bases.")

        else:

            length = len(seq)

            if seq_type == "dsDNA":
                mw_per_bp = 650
            else:
                mw_per_bp = 330

            mol_weight = length * mw_per_bp

            st.subheader("Sequence Properties")

            st.write(f"Length: **{length} bp**")
            st.write(f"Molecular weight: **{mol_weight:.2e} g/mol**")

            pmol_per_ng = (1e-9 / mol_weight) * 1e12

            st.write(f"1 ng corresponds to **{pmol_per_ng:.4f} pmol**")

            if mass_ng > 0:

                moles = (mass_ng * 1e-9) / mol_weight
                pmol = moles * 1e12

                st.subheader("Mass Conversion")

                st.write(f"{mass_ng} ng = **{pmol:.4f} pmol**")


# ============================================
# 9. SEQUENCE TOOLS
# ============================================

elif st.session_state.page == "tools":

    st.markdown("### Sequence Analysis Tools")

    # ---------------------------
    # COMPLEMENT
    # ---------------------------

    seq = st.text_area("Paste DNA sequence")

    if seq:

        seq = seq.upper().replace("\n","").replace(" ","")

        comp_map = str.maketrans("ATGC","TACG")

        complement = seq.translate(comp_map)
        rev = complement[::-1]

        st.subheader("Complement")
        st.code(complement)

        st.subheader("Reverse Complement")
        st.code(rev)

    st.markdown("---")

# ---------------------------
# ALIGNMENT (ADVANCED)
# ---------------------------
    # ---------------------------
# ALIGNMENT (ADVANCED)
# ---------------------------

st.subheader("Pairwise Alignment")

mode = st.selectbox(
    "Alignment type",
    ["Local (Smith-Waterman)", "Global (Affine gaps)"],
    key="align_mode"
)

s1 = st.text_area("First sequence", key="align_s1")
s2 = st.text_area("Second sequence", key="align_s2")

# ---------------------------
# SCORING (TRANSITION AWARE)
# ---------------------------

transitions = {"AG", "GA", "CT", "TC"}

def score_pair(a, b):
    if a == b:
        return 2
    elif a + b in transitions:
        return -0.5   # transition (less severe)
    else:
        return -1     # transversion


# ---------------------------
# VISUALIZATION (COLORED)
# ---------------------------

def render_alignment_colored(a1, a2):
    html = ""
    for x, y in zip(a1, a2):
        if x == y:
            html += f"<span style='color:green'>{x}</span>"
        elif "-" in (x, y):
            html += f"<span style='color:gray'>{x}</span>"
        else:
            html += f"<span style='color:red'>{x}</span>"

    html2 = ""
    for x, y in zip(a1, a2):
        if x == y:
            html2 += f"<span style='color:green'>{y}</span>"
        elif "-" in (x, y):
            html2 += f"<span style='color:gray'>{y}</span>"
        else:
            html2 += f"<span style='color:red'>{y}</span>"

    st.markdown(html, unsafe_allow_html=True)
    st.markdown(html2, unsafe_allow_html=True)


def render_alignment_plain(a1, a2):
    visual = "".join(
        "|" if x == y else " " if "-" in (x, y) else "."
        for x, y in zip(a1, a2)
    )
    st.code(f"{a1}\n{visual}\n{a2}")


# ---------------------------
# SMITH-WATERMAN
# ---------------------------

def smith_waterman(s1, s2):

    gap = -2
    n, m = len(s1), len(s2)

    score = [[0]*(m+1) for _ in range(n+1)]

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, n+1):
        for j in range(1, m+1):

            diag = score[i-1][j-1] + score_pair(s1[i-1], s2[j-1])
            delete = score[i-1][j] + gap
            insert = score[i][j-1] + gap

            score[i][j] = max(0, diag, delete, insert)

            if score[i][j] > max_score:
                max_score = score[i][j]
                max_pos = (i, j)

    align1, align2 = "", ""
    i, j = max_pos
    matches, length = 0, 0

    while i > 0 and j > 0 and score[i][j] != 0:

        if score[i][j] == score[i-1][j-1] + score_pair(s1[i-1], s2[j-1]):
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            if s1[i-1] == s2[j-1]:
                matches += 1
            i -= 1
            j -= 1

        elif score[i][j] == score[i-1][j] - 2:
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        else:
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            j -= 1

        length += 1

    identity = (matches / length * 100) if length > 0 else 0

    return align1, align2, max_score, identity


# ---------------------------
# AFFINE GAP (GOTOH)
# ---------------------------

def affine_global(s1, s2):

    gap_open, gap_extend = -5, -1

    n, m = len(s1), len(s2)

    M = [[0]*(m+1) for _ in range(n+1)]
    X = [[0]*(m+1) for _ in range(n+1)]
    Y = [[0]*(m+1) for _ in range(n+1)]

    for i in range(1, n+1):
        M[i][0] = gap_open + (i-1)*gap_extend
        X[i][0] = M[i][0]

    for j in range(1, m+1):
        M[0][j] = gap_open + (j-1)*gap_extend
        Y[0][j] = M[0][j]

    for i in range(1, n+1):
        for j in range(1, m+1):

            X[i][j] = max(
                M[i-1][j] + gap_open,
                X[i-1][j] + gap_extend
            )

            Y[i][j] = max(
                M[i][j-1] + gap_open,
                Y[i][j-1] + gap_extend
            )

            M[i][j] = max(
                M[i-1][j-1] + score_pair(s1[i-1], s2[j-1]),
                X[i][j],
                Y[i][j]
            )

    align1, align2 = "", ""
    i, j = n, m

    while i > 0 or j > 0:

        if i > 0 and j > 0 and M[i][j] == M[i-1][j-1] + score_pair(s1[i-1], s2[j-1]):
            align1 = s1[i-1] + align1
            align2 = s2[j-1] + align2
            i -= 1
            j -= 1

        elif i > 0 and M[i][j] == X[i][j]:
            align1 = s1[i-1] + align1
            align2 = "-" + align2
            i -= 1

        else:
            align1 = "-" + align1
            align2 = s2[j-1] + align2
            j -= 1

    return align1, align2, M[n][m]


# ---------------------------
# EXPORT FUNCTIONS
# ---------------------------

def export_fasta(a1, a2):
    return f">Sequence_1\n{a1}\n>Sequence_2\n{a2}"


def export_clustal(a1, a2):
    return f"CLUSTAL format\n\nSeq1    {a1}\nSeq2    {a2}"


# ---------------------------
# RUN ALIGNMENT
# ---------------------------

if st.button("Run Alignment"):

    s1 = s1.upper().replace("\n","").replace(" ","")
    s2 = s2.upper().replace("\n","").replace(" ","")

    valid_bases = set("ATGC")

    if not s1 or not s2:
        st.warning("Please provide both sequences.")
        st.stop()

    if not set(s1).issubset(valid_bases) or not set(s2).issubset(valid_bases):
        st.error("Sequences must contain only A, T, G, C.")
        st.stop()

    if mode == "Local (Smith-Waterman)":
        a1, a2, score, identity = smith_waterman(s1, s2)
        st.write(f"Max Alignment Score: {score}")
        st.write(f"Identity: {identity:.2f}%")

    else:
        a1, a2, score = affine_global(s1, s2)
        st.write(f"Alignment Score: {score}")

        mismatch = any(a != b for a, b in zip(a1[-5:], a2[-5:]))
        if mismatch:
            st.warning("⚠ Possible 3′ primer mismatch detected")
        else:
            st.success("3′ end appears well matched")

    # Render
    render_alignment_plain(a1, a2)
    st.markdown("### Colored Alignment")
    render_alignment_colored(a1, a2)

    # Export
    st.markdown("### Export Alignment")

    fasta = export_fasta(a1, a2)
    clustal = export_clustal(a1, a2)

    st.download_button("Download FASTA", fasta, file_name="alignment.fasta")
    st.download_button("Download CLUSTAL", clustal, file_name="alignment.aln")


    st.markdown("---")
# ---------------------------
# PRIMER GENERATOR
# ---------------------------
    st.subheader("Primer Generator")
    template = st.text_area("Template DNA")
    length = st.number_input("Primer length", value=20)

    if st.button("Generate Primers"):

        template = template.upper().replace("\n","").replace(" ","")

        comp = str.maketrans("ATGC","TACG")

        forward = template[:length]
        reverse = template[-length:].translate(comp)[::-1]

        def gc(seq):
            return (seq.count("G")+seq.count("C"))/len(seq)*100

        def tm(seq):
            return 2*(seq.count("A")+seq.count("T"))+4*(seq.count("G")+seq.count("C"))

        def gc_clamp(seq):
            last5=seq[-5:]
            gc_count=last5.count("G")+last5.count("C")
            if gc_count>=2:
                return "Good GC clamp"
            return "Weak GC clamp"

        def reverse_comp(seq):
            return seq.translate(comp)[::-1]

        def hairpin_check(seq):
            rc=reverse_comp(seq)
            for i in range(len(seq)-4):
                if seq[i:i+4] in rc:
                    return "Possible hairpin"
            return "No strong hairpin detected"

        def dimer_check(p1,p2):
            rc2=reverse_comp(p2)
            for i in range(len(p1)-4):
                if p1[i:i+4] in rc2:
                    return "Possible primer dimer"
            return "No strong dimer detected"

        st.write("Forward Primer")
        st.code(forward)
        st.write(f"GC%: {gc(forward):.1f}")
        st.write(f"Tm: {tm(forward)} °C")
        st.write(gc_clamp(forward))
        st.write(hairpin_check(forward))

        st.write("Reverse Primer")
        st.code(reverse)
        st.write(f"GC%: {gc(reverse):.1f}")
        st.write(f"Tm: {tm(reverse)} °C")
        st.write(gc_clamp(reverse))
        st.write(hairpin_check(reverse))

        st.write("Dimer Check")
        st.write(dimer_check(forward,reverse))


# ============================================
# 10. FOOTER
# ============================================

st.markdown("---")
st.caption("🧬 RNMolBioTools • Molecular Biology Utilities")
