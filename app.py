# --------------------------------------------
# 1. IMPORT LIBRARIES
# --------------------------------------------

import streamlit as st


# --------------------------------------------
# 2. PAGE CONFIGURATION
# --------------------------------------------

st.set_page_config(
    page_title="DNA Cloning Calculator | Njue BioTools",
    page_icon="🧬",
)

st.markdown(
"""
<div style="
background: linear-gradient(90deg,#2C3E50,#4CA1AF);
padding:14px;
border-radius:8px;
display:flex;
align-items:center;
gap:12px;
color:white;
box-shadow:0px 2px 6px rgba(0,0,0,0.2);
">

<div style="font-size:32px;">
🧬
</div>

<div>
<div style="font-size:18px;font-weight:600;">
DNA Cloning Calculator
</div>

<div style="font-size:13px;opacity:0.85;">
Njue BioTools
</div>
</div>

</div>
""",
unsafe_allow_html=True
)


# --------------------------------------------
# 3. HEADER (Njue BioTools styled background)
# --------------------------------------------

st.markdown("your styled header code here", unsafe_allow_html=True)


# --------------------------------------------
# 4. NAVIGATION BUTTONS  ← PLACE CODE HERE
# --------------------------------------------

if "page" not in st.session_state:
    st.session_state.page = "cloning"

col1, col2, col3 = st.columns(3)

with col1:
    if st.button("🧬 DNA Cloning Calculator", use_container_width=True):
        st.session_state.page = "cloning"

with col2:
    if st.button("🔢 DNA Copy Number", use_container_width=True):
        st.session_state.page = "copy"

with col3:
    if st.button("⚖ DNA/RNA Mass Converter", use_container_width=True):
        st.session_state.page = "conversion"

st.markdown("---")


# --------------------------------------------
# 5. SHOW SELECTED CALCULATOR
# --------------------------------------------

if st.session_state.page == "cloning":

    st.markdown("#### DNA Cloning Calculator")

import streamlit as st

# --------------------------------------------
# DNA Cloning Calculator
# --------------------------------------------

st.markdown(
"""
<h4 style='margin-bottom:5px;'>DNA Cloning Calculator</h4>
""",
unsafe_allow_html=True
)

vector_kb = st.number_input("Vector DNA size (kb)", min_value=0.0)
insert_kb = st.number_input("Insert DNA size (kb)", min_value=0.0)
ratio = st.text_input("Vector : Insert ratio; You can adjust", "1:4")
total_ng = st.number_input("Total DNA (in ng) e.g. 100 ng, user adjust accordingly", value=100.0)

if st.button("Calculate"):

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

# --------------------------------------------
# DNA Cloning Calculator info
# --------------------------------------------

st.caption("ⓘ Calculation assumptions")

with st.expander("ℹ️ How the ligation calculation works"):

    st.markdown("""
**Molar ratios**

DNA ligation efficiency depends on the **molar ratio of vector to insert**, not their mass.  
The calculator converts DNA mass (ng) into **moles** to maintain the correct ratio.

---

**Average molecular weight of dsDNA**

1 base pair (bp) of double-stranded DNA has an average molecular weight of:

**≈ 650–660 g/mol per bp**

This calculator uses **650 g/mol per bp**.

---

**Converting DNA mass to pmol**

The conversion used is:

pmol = (ng × 1000) / (bp × 650)

Where:

- **ng × 1000** converts nanograms → picograms
- **bp × 650** gives the molecular weight of the DNA fragment

---

**Vector and insert mass distribution**

The total DNA amount is distributed between vector and insert according to:

vector_ng = total_ng / (1 + ratio × insert_bp / vector_bp)

The insert mass is:

insert_ng = total_ng − vector_ng

---

**Typical ligation ratios**

Common molar ratios used in cloning:

- **1 : 3** → sticky-end ligations
- **1 : 4–6** → blunt-end ligations
- **1 : 1** → large fragment ligations
""")
    
# --------------------------------------------
# dsDNA Copy Number Calculator 
# --------------------------------------------

elif st.session_state.page == "copy":

    st.markdown("#### dsDNA Copy Number Calculator")

st.markdown(
"""
<h4 style='margin-bottom:5px;'>dsDNA Copy Number Calculator</h4>
""",
unsafe_allow_html=True
)

dna_kb = st.number_input("DNA length (kb)", min_value=0.0, key="dna_kb")
dna_ng = st.number_input("DNA amount (ng)", min_value=0.0, key="dna_ng")

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
    
# --------------------------------------------
# DNA Copy number calculation info
# --------------------------------------------

st.caption("ⓘ Calculation assumptions")

with st.expander("ℹ️ How DNA copy number is calculated"):

    st.markdown("""
**Average molecular weight of dsDNA**

1 base pair (bp) of double-stranded DNA has an average molecular weight of:

**≈ 650–660 g/mol**

This calculator assumes **650 g/mol per base pair**.

---

**Molar mass of a DNA molecule**

Molar mass is calculated as:

Molar mass = bp × 650 g/mol

Example:

A **4 kb plasmid**

= 4000 bp × 650  
= **2.6 × 10⁶ g/mol**

---

**Converting DNA mass to moles**

DNA mass is first converted from **ng to grams**:

moles = (DNA mass in g) / molar mass

---

**Calculating number of DNA molecules**

The number of molecules is determined using **Avogadro's constant**:

6.022 × 10²³ molecules/mol

copies = moles × 6.022 × 10²³

---

**Example**

50 ng of a 4 kb plasmid:

Molar mass = 2.6 × 10⁶ g/mol  
Moles ≈ 1.9 × 10⁻¹⁴ mol  

DNA copies ≈ **1.1 × 10¹⁰ molecules**
""")

# --------------------------------------------
# DNA/RNA Mass Conversion Calculator
# --------------------------------------------

elif st.session_state.page == "conversion":

    st.markdown("#### DNA/RNA Mass Conversion Tool")

st.markdown("---")
st.markdown("#### DNA/RNA Mass Unit Converter")

na_type = st.selectbox(
    "Molecule type",
    ["dsDNA", "ssDNA", "RNA"]
)

col1, col2, col3, col4, col5 = st.columns(5)

g = col1.number_input("g", value=0.0)
mg = col2.number_input("mg", value=0.0)
ug = col3.number_input("µg", value=0.0)
ng = col4.number_input("ng", value=0.0)
pg = col5.number_input("pg", value=0.0)

# determine which field user changed
if g > 0:
    mg = g * 1e3
    ug = g * 1e6
    ng = g * 1e9
    pg = g * 1e12

elif mg > 0:
    g = mg / 1e3
    ug = mg * 1e3
    ng = mg * 1e6
    pg = mg * 1e9

elif ug > 0:
    g = ug / 1e6
    mg = ug / 1e3
    ng = ug * 1e3
    pg = ug * 1e6

elif ng > 0:
    g = ng / 1e9
    mg = ng / 1e6
    ug = ng / 1e3
    pg = ng * 1e3

elif pg > 0:
    g = pg / 1e12
    mg = pg / 1e9
    ug = pg / 1e6
    ng = pg / 1e3

st.write("Converted values updated automatically.")
    
st.markdown("---")
st.caption("🧬 Njue BioTools • Molecular Biology Utilities")
