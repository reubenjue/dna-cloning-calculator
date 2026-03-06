import streamlit as st

col1, col2 = st.columns([1,6])

with col1:
    st.markdown("## 🧬")

with col2:
    st.markdown("### DNA Cloning Calculator")
    st.caption("Njue BioTools")
    
import streamlit as st

st.title("DNA Cloning Calculator")

vector_kb = st.number_input("Vector DNA size (kb)", min_value=0.0)
insert_kb = st.number_input("Insert DNA size (kb)", min_value=0.0)
ratio = st.text_input("Vector : Insert ratio; You can adjust", "1:4")
total_ng = st.number_input("Total DNA (in ng); 100 ng is recommended for ligation, but user can adjust accordingly", value=100.0)

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
    
st.markdown("---")
st.header("dsDNA Copy Number Calculator")

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
