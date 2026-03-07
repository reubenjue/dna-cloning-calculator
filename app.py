# ============================================
# 1. IMPORT LIBRARIES
# ============================================

import streamlit as st
import matplotlib.pyplot as plt
import numpy as np


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
# 6. COPY NUMBER CALCULATOR
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

        st.write(f"Molar Mass: {molar_mass:.2e}")
        st.write(f"Moles: {moles:.2e}")
        st.write(f"pmol: {pmol:.4f}")
        st.write(f"Copies: {copies:.2e}")

    st.caption("ⓘ Calculation assumptions")

    with st.expander("How DNA copy number is calculated"):

        st.markdown("""
Average molecular weight:

**650 g/mol per base pair**

Copies calculated using **Avogadro's number**:

6.022 × 10²³ molecules/mol
""")


# ============================================
# 7. MASS CONVERSION
# ============================================

elif st.session_state.page == "conversion":

    st.markdown("### DNA/RNA Mass Unit Converter")

    g = st.number_input("g")
    st.write(f"mg: {g*1e3}")
    st.write(f"µg: {g*1e6}")
    st.write(f"ng: {g*1e9}")
    st.write(f"pg: {g*1e12}")


# ============================================
# 8. SEQUENCE MASS
# ============================================

elif st.session_state.page == "sequence":

    st.markdown("### DNA Sequence Mass Calculator")

    seq = st.text_area("Paste DNA sequence")

    if seq:

        seq = seq.upper().replace("\n","")

        length=len(seq)
        mol_weight=length*650

        st.write(f"Length: {length}")
        st.write(f"Molecular weight: {mol_weight:.2e}")


# ============================================
# 9. SEQUENCE TOOLS
# ============================================

elif st.session_state.page == "tools":

    st.markdown("### Sequence Tools")

    template = st.text_area("Template DNA")

    if template:

        template=template.upper().replace("\n","")

        st.subheader("Restriction Enzyme Site Scanner")

        enzymes={
            "EcoRI":"GAATTC",
            "BamHI":"GGATCC",
            "HindIII":"AAGCTT",
            "XhoI":"CTCGAG"
        }

        cut_positions=[]

        for enzyme,site in enzymes.items():

            start=0
            while True:
                pos=template.find(site,start)
                if pos==-1:
                    break
                st.write(f"{enzyme} ({site}) at position {pos+1}")
                cut_positions.append((enzyme,pos+1))
                start=pos+1


        # ----------------------------------------
        # Restriction Digest Simulator
        # ----------------------------------------

        st.subheader("Restriction Digest Simulator")

        if cut_positions:

            cuts=sorted([p for _,p in cut_positions])
            fragments=[]
            prev=0

            for c in cuts:
                fragments.append(c-prev)
                prev=c

            fragments.append(len(template)-prev)

            st.write("Predicted fragment sizes (bp):")
            st.write(sorted(fragments))


        # ----------------------------------------
        # Graphical Restriction Map
        # ----------------------------------------

        st.subheader("Plasmid Restriction Map")

        seq_len=len(template)

        fig, ax = plt.subplots(figsize=(6,6))

        circle = plt.Circle((0,0),1,fill=False,linewidth=2)
        ax.add_patch(circle)

        for enzyme,pos in cut_positions:

            angle = 2*np.pi*(pos/seq_len)

            x = np.cos(angle)
            y = np.sin(angle)

            ax.plot([0,x],[0,y])

            ax.text(1.2*x,1.2*y,f"{enzyme}\n{pos}",ha='center')

        ax.set_xlim(-1.5,1.5)
        ax.set_ylim(-1.5,1.5)
        ax.axis("off")

        st.pyplot(fig)


# ============================================
# FOOTER
# ============================================

st.markdown("---")
st.caption("🧬 RNMolBioTools • Molecular Biology Utilities")
