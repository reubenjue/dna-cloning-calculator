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
# 5. DNA CLONING CALCULATOR (UNCHANGED)
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

# ============================================
# 6. DNA COPY NUMBER CALCULATOR (UNCHANGED)
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


# ============================================
# 7. DNA / RNA MASS CONVERSION (UNCHANGED)
# ============================================

elif st.session_state.page == "conversion":

    st.markdown("### DNA/RNA Mass Unit Converter")

    na_type = st.selectbox(
        "Molecule type",
        ["dsDNA", "ssDNA", "RNA"]
    )

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


# ============================================
# 8. DNA SEQUENCE MASS CALCULATOR (UNCHANGED)
# ============================================

elif st.session_state.page == "sequence":

    st.markdown("### DNA Sequence Mass Calculator")

    seq_type = st.selectbox(
        "Sequence Type",
        ["dsDNA", "ssDNA"]
    )

    sequence = st.text_area(
        "Paste DNA sequence (A, T, G, C only)",
        height=150
    )

    mass_ng = st.number_input(
        "Optional: DNA amount (ng)",
        min_value=0.0,
        value=0.0
    )

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
# 9. SEQUENCE TOOLS (NEW FEATURES)
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
    # ALIGNMENT
    # ---------------------------

    st.subheader("Pairwise Alignment")

    s1 = st.text_area("First sequence")
    s2 = st.text_area("Second sequence")

    if st.button("Run Alignment"):

        s1 = s1.upper().replace("\n","").replace(" ","")
        s2 = s2.upper().replace("\n","").replace(" ","")

        match = 1
        mismatch = -1
        gap = -1

        n=len(s1)
        m=len(s2)

        score=[[0]*(m+1) for _ in range(n+1)]

        for i in range(n+1):
            score[i][0]=i*gap

        for j in range(m+1):
            score[0][j]=j*gap

        for i in range(1,n+1):
            for j in range(1,m+1):

                diag=score[i-1][j-1]+(match if s1[i-1]==s2[j-1] else mismatch)
                delete=score[i-1][j]+gap
                insert=score[i][j-1]+gap

                score[i][j]=max(diag,delete,insert)

        align1=""
        align2=""

        i=n
        j=m
        matches=0
        length=0

        while i>0 or j>0:

            if i>0 and j>0 and score[i][j]==score[i-1][j-1]+(match if s1[i-1]==s2[j-1] else mismatch):

                align1=s1[i-1]+align1
                align2=s2[j-1]+align2

                if s1[i-1]==s2[j-1]:
                    matches+=1

                i-=1
                j-=1

            elif i>0 and score[i][j]==score[i-1][j]+gap:

                align1=s1[i-1]+align1
                align2="-"+align2
                i-=1

            else:

                align1="-"+align1
                align2=s2[j-1]+align2
                j-=1

            length+=1

        percent=(matches/length)*100

        visual=""

        for a,b in zip(align1,align2):

            if a==b:
                visual+="|"
            elif a=="-" or b=="-":
                visual+=" "
            else:
                visual+="."

        mismatch=False

        for a,b in zip(align1[-5:],align2[-5:]):
            if a!=b:
                mismatch=True

        if mismatch:
            warn="⚠ Possible 3′ primer mismatch detected"
        else:
            warn="3′ end appears well matched"

        st.write(f"Alignment Score: {score[n][m]}")
        st.write(f"Complementarity: {percent:.2f}%")
        st.write(warn)

        st.code(f"{align1}\n{visual}\n{align2}")

    st.markdown("---")

    # ---------------------------
    # PRIMER DESIGN
    # ---------------------------

    st.subheader("Primer Generator")

    template = st.text_area("Template DNA")
    length = st.number_input("Primer length", value=20)

    if st.button("Generate Primers"):

        template = template.upper().replace("\n","").replace(" ","")

        forward = template[:length]

        comp = str.maketrans("ATGC","TACG")

        reverse = template[-length:].translate(comp)[::-1]

        def gc(seq):
            return (seq.count("G")+seq.count("C"))/len(seq)*100

        def tm(seq):
            return 2*(seq.count("A")+seq.count("T"))+4*(seq.count("G")+seq.count("C"))

        st.write("Forward Primer")
        st.code(forward)
        st.write(f"GC%: {gc(forward):.1f}")
        st.write(f"Tm: {tm(forward)} °C")

        st.write("Reverse Primer")
        st.code(reverse)
        st.write(f"GC%: {gc(reverse):.1f}")
        st.write(f"Tm: {tm(reverse)} °C")


# ============================================
# FOOTER
# ============================================

st.markdown("---")
st.caption("🧬 RNMolBioTools • Molecular Biology Utilities")
