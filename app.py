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
