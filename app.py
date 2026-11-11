# ============================================================================
# RNMolBioTools — DNA Cloning Calculator, V9 (Streamlit)
# All scientific logic lives in molbio_core.py (same folder).
# ============================================================================

import pandas as pd
import streamlit as st

import molbio_core as mb

st.set_page_config(page_title="DNA Cloning Calculator | RNMolBioTools",
                   page_icon="🧬", layout="wide")

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------


def sci(x, digits=3):
    """Readable number: plain for 0.01–99,999, scientific otherwise."""
    if x is None:
        return "–"
    if x == 0:
        return "0"
    ax = abs(x)
    if 0.01 <= ax < 1e5:
        return f"{x:,.{digits}g}" if ax < 1 else f"{x:,.{max(digits - 1, 1)}f}".rstrip("0").rstrip(".")
    return f"{x:.{digits}e}"


def show_problems(problems):
    for p in problems:
        st.error(p)
    return bool(problems)


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


MOLECULES = ["dsDNA", "ssDNA", "RNA"]
LEN_UNIT = {"dsDNA": "bp", "ssDNA": "nt", "RNA": "nt"}

# ----------------------------------------------------------------------------
# Header
# ----------------------------------------------------------------------------

st.markdown(
    """
<div style="background:linear-gradient(90deg,#2C3E50,#4CA1AF);padding:16px;border-radius:8px;
display:flex;align-items:center;gap:12px;color:white;box-shadow:0 2px 6px rgba(0,0,0,.2);">
<div style="font-size:34px;">🧬</div>
<div><div style="font-size:20px;font-weight:600;">DNA Cloning Calculator</div>
<div style="font-size:13px;opacity:.9;">RNMolBioTools · V9</div></div></div>
""",
    unsafe_allow_html=True,
)
st.write("")

# ----------------------------------------------------------------------------
# Sidebar: global conventions and PCR conditions
# ----------------------------------------------------------------------------

with st.sidebar:
    st.header("Settings")
    mw_method = st.radio(
        "Molecular-weight convention", mb.MW_METHODS, index=0,
        help="Exact: dsDNA = 617.96·bp + 36.04; ssDNA = 308.97·nt + 18.02; "
             "RNA = 321.47·nt + 18.02 g/mol (NEB/QIAGEN). Legacy: 650 / 330 / 340 g/mol "
             "per bp or nt (sodium-salt rule of thumb, ~5 % higher).")
    st.subheader("PCR conditions")
    st.caption("Used for Tm and ΔG in the primer tools.")
    na = st.number_input("Monovalent cations Na⁺/K⁺ (mM)", 0.0, 1000.0, 50.0, 5.0)
    mg = st.number_input("Mg²⁺ (mM)", 0.0, 20.0, 1.5, 0.25)
    dntp = st.number_input("Total dNTP (mM)", 0.0, 10.0, 0.8, 0.1,
                           help="Sum of all four dNTPs (0.2 mM each = 0.8 mM). dNTPs chelate Mg²⁺.")
    primer_nm = st.number_input("Primer concentration (nM, each)", 1.0, 5000.0, 250.0, 50.0)
    cond = mb.PCRConditions(na, mg, dntp, primer_nm)
    st.caption(f"Na⁺-equivalent: {cond.na_equivalent_M() * 1000:.0f} mM · "
               f"ΔG engine: {mb.thermo_backend()}")

tabs = st.tabs(["🧬 Ligation & Assembly", "🔢 Copy Number & Standards",
                "⚖ Mass ↔ Moles", "🧾 Sequence Properties",
                "🧪 Sequence Tools", "🎯 Primers"])

# ============================================================================
# 1. LIGATION & ASSEMBLY
# ============================================================================

with tabs[0]:
    st.markdown("### Ligation & Assembly Calculator")
    rxn = st.radio("Reaction type",
                   ["T4 DNA ligation (sticky / blunt ends)",
                    "Seamless assembly (Gibson / NEBuilder HiFi / In-Fusion)"],
                   horizontal=True)
    is_asm = rxn.startswith("Seamless")

    c1, c2 = st.columns(2)
    with c1:
        vector_bp = st.number_input("Linearised vector length (bp)", min_value=1, value=4000, step=100)
        n_ins = st.number_input("Number of inserts", min_value=1, max_value=6 if is_asm else 3,
                                value=1, step=1)
        inserts = []
        for k in range(int(n_ins)):
            a, b = st.columns(2)
            bp = a.number_input(f"Insert {k + 1} length (bp)", min_value=1, value=1000,
                                step=50, key=f"ins_bp_{k}")
            r = b.number_input(f"Insert {k + 1} : vector molar ratio", min_value=0.1,
                               value=2.0 if is_asm else 3.0, step=0.5, key=f"ins_r_{k}",
                               help="Moles of this insert per mole of vector (3 = 3-fold excess).")
            inserts.append((int(bp), float(r)))
    with c2:
        mode_labels = {"Fix vector mass (ng) — recommended": "vector_ng",
                       "Fix total DNA mass (ng)": "total_ng",
                       "Fix total DNA amount (pmol)": "total_pmol"}
        default_mode = 2 if is_asm else 0
        mode_lbl = st.selectbox("How to set the amount", list(mode_labels), index=default_mode)
        mode = mode_labels[mode_lbl]
        default_amt = {"vector_ng": 50.0, "total_ng": 100.0, "total_pmol": 0.1}[mode]
        amount = st.number_input(
            {"vector_ng": "Vector (ng)", "total_ng": "Total DNA (ng)",
             "total_pmol": "Total DNA (pmol)"}[mode],
            min_value=0.0001, value=default_amt, step=0.01 if mode == "total_pmol" else 5.0,
            format="%.4f" if mode == "total_pmol" else "%.1f", key=f"amt_{mode}")
        rxn_ul = st.number_input("Reaction volume (µL)", min_value=1.0, value=20.0, step=1.0)
        other_ul = st.number_input(
            "Volume of other components (µL)", min_value=0.0,
            value=10.0 if is_asm else 3.0, step=0.5,
            help="Buffer + enzyme. T4: 2 µL 10× buffer + 1 µL ligase. "
                 "HiFi/Gibson 2× master mix: half the reaction volume.")
        use_stocks = st.checkbox("I know my stock concentrations → give pipetting volumes")
        stocks = None
        if use_stocks:
            stocks = [st.number_input("Vector stock (ng/µL)", min_value=0.001, value=50.0,
                                      key="stock_v")]
            for k in range(int(n_ins)):
                stocks.append(st.number_input(f"Insert {k + 1} stock (ng/µL)", min_value=0.001,
                                              value=20.0, key=f"stock_{k}"))

    try:
        res = mb.ligation_mix(int(vector_bp), inserts, mode, float(amount), stocks=stocks,
                              reaction_ul=rxn_ul, other_ul=other_ul)
        rows = []
        for c in res.components:
            rows.append({"Component": c.name, "Length (bp)": c.bp,
                         "Molar parts": c.ratio, "Mass (ng)": round(c.ng, 2),
                         "Amount (fmol)": round(c.pmol * 1000, 2),
                         "Volume (µL)": None if c.ul is None else round(c.ul, 2)})
        st.subheader("Reaction mix")
        st.dataframe(pd.DataFrame(rows), hide_index=True)
        m1, m2, m3 = st.columns(3)
        m1.metric("Total DNA", f"{res.total_ng:.1f} ng")
        m2.metric("Total DNA", f"{res.total_pmol:.4f} pmol")
        if res.dna_volume_ul is not None:
            water = rxn_ul - other_ul - res.dna_volume_ul
            m3.metric("Water to add", f"{max(water, 0):.2f} µL")
        for w in res.warnings:
            st.warning(w)
        if is_asm:
            n_frag = 1 + int(n_ins)
            lo, hi = (0.03, 0.2) if n_frag <= 3 else (0.2, 0.5)
            if not lo <= res.total_pmol <= hi:
                st.info(f"For {n_frag} fragments, NEBuilder HiFi recommends "
                        f"{lo}–{hi} pmol total DNA (check your kit's manual).")
        else:
            if not 0.01 <= res.components[0].pmol <= 0.1:
                st.info("Typical T4 ligations use about 0.01–0.1 pmol vector "
                        "(e.g. 50 ng of a 4 kb vector ≈ 0.02 pmol).")
    except ValueError as e:
        st.error(str(e))

    with st.expander("Calculation and guidance"):
        st.markdown(r"""
**Moles from mass:** pmol = ng × 1000 / MW, with MW(dsDNA, linear) = 617.96·bp + 36.04 g/mol.

**Fix vector mass (recommended):** insert ng = ratio × vector pmol × MW(insert) / 1000.
This keeps the vector — and therefore the number of clonable ends — constant, which is
what controls ligation efficiency.

**Fix total mass:** vector ng = total / (1 + Σ ratioᵢ · MWᵢ / MW(vector)).

**Starting points (always optimise empirically):**
- T4 ligation: 1:3 vector:insert; try 1:1 to 1:10. Very short inserts (< ~200 bp)
  usually need higher excess; blunt-end ligations are less efficient — increase ligase
  or add PEG rather than only raising the ratio.
- NEBuilder HiFi / Gibson: 2–3 fragments 0.03–0.2 pmol total, vector:insert 1:2;
  4–6 fragments 0.2–0.5 pmol, 1:1; inserts < 200 bp 1:5.
- Keep total DNA modest (≈ 1–10 ng/µL); high concentrations favour concatemers.
""")

# ============================================================================
# 2. COPY NUMBER & qPCR STANDARDS
# ============================================================================

with tabs[1]:
    st.markdown("### Copy Number Calculator")
    kind = st.radio("Template", ["Plasmid / PCR product / transcript", "Genomic DNA"],
                    horizontal=True)
    if kind.startswith("Plasmid"):
        c1, c2 = st.columns(2)
        with c1:
            mol = st.selectbox("Molecule", MOLECULES, key="cn_mol")
            length = st.number_input(f"Length ({LEN_UNIT[mol]})", min_value=1, value=4000,
                                     step=100, key="cn_len")
            circ = mol == "dsDNA" and st.checkbox("Circular (uncut plasmid)", key="cn_circ")
        with c2:
            how = st.radio("Input", ["Concentration (ng/µL)", "Mass (ng)"], key="cn_how")
            val = st.number_input(how, min_value=0.0, value=10.0, key="cn_val")
        try:
            mw = mb.mw_from_length(int(length), mol, circ, mw_method)
            m1, m2, m3 = st.columns(3)
            m1.metric("Molecular weight", f"{sci(mw)} g/mol")
            if how.startswith("Concentration"):
                m2.metric("Copies / µL", sci(mb.copies_from_ng(val, mw)))
                m3.metric("Molarity", f"{sci(mb.ng_per_ul_to_nM(val, mw))} nM")
                stock_copies = mb.copies_from_ng(val, mw)
            else:
                m2.metric("Copies", sci(mb.copies_from_ng(val, mw)))
                m3.metric("Amount", f"{sci(mb.ng_to_pmol(val, mw))} pmol")
                stock_copies = None
        except ValueError as e:
            st.error(str(e))
            stock_copies = None
    else:
        c1, c2 = st.columns(2)
        genome_mb = c1.number_input("Haploid genome size (Mb)", min_value=0.001, value=3100.0,
                                    help="Human ≈ 3,100 Mb; mouse ≈ 2,700 Mb; E. coli ≈ 4.6 Mb.")
        gdna = c2.number_input("Genomic DNA (ng)", min_value=0.0, value=10.0)
        pg_per = mb.haploid_genome_mass_pg(genome_mb * 1e6)
        m1, m2 = st.columns(2)
        m1.metric("Mass per haploid genome", f"{sci(pg_per)} pg")
        m2.metric("Haploid genome copies", sci(gdna * 1000 / pg_per))
        st.caption("A single-copy gene is present once per haploid genome; a diploid cell "
                   "contains two haploid genome equivalents.")
        stock_copies = None

    st.markdown("---")
    st.markdown("#### Serial dilution planner (qPCR standard curve)")
    d1, d2, d3, d4, d5 = st.columns(5)
    stock = d1.number_input("Stock (copies/µL)", min_value=1.0,
                            value=float(stock_copies) if stock_copies else 1e10, format="%.3e")
    top = d2.number_input("Top standard (copies/µL)", min_value=1e-3, value=1e7, format="%.3e")
    fold = d3.number_input("Fold dilution", min_value=1.5, value=10.0, step=1.0)
    npts = d4.number_input("Points", min_value=1, max_value=12, value=6, step=1)
    vol = d5.number_input("Volume kept per tube (µL)", min_value=1.0, value=50.0)
    try:
        rows, warns = mb.serial_dilution(stock, top, fold, int(npts), vol)
        df = pd.DataFrame(rows)
        df["Concentration"] = df["Concentration"].map(lambda x: f"{x:.2e} copies/µL")
        for col in ("Transfer (µL)", "Diluent (µL)", "Total made (µL)"):
            df[col] = df[col].map(lambda x: round(x, 3))
        st.dataframe(df, hide_index=True)
        for w in warns:
            st.warning(w)
        st.caption("Use carrier (e.g. 10 ng/µL yeast tRNA or salmon sperm DNA) in the diluent "
                   "for low-copy standards to limit adsorption losses.")
    except ValueError as e:
        st.error(str(e))

# ============================================================================
# 3. MASS ↔ MOLES CONVERTER
# ============================================================================

with tabs[2]:
    st.markdown("### Mass ↔ Moles ↔ Molarity")
    st.caption("Mass-unit conversion alone does not depend on the molecule; converting "
               "to moles does — that is what molecule type and length are used for.")
    c1, c2, c3 = st.columns(3)
    with c1:
        mol = st.selectbox("Molecule", MOLECULES, key="cv_mol")
        length = st.number_input(f"Length ({LEN_UNIT[mol]})", min_value=1, value=1000,
                                 step=100, key="cv_len")
    with c2:
        qty = st.number_input("Value", min_value=0.0, value=100.0, format="%.6g", key="cv_val")
        unit = st.selectbox("Unit", list(mb.MASS_UNITS) + list(mb.MOLE_UNITS), index=3)
    with c3:
        conc = st.number_input("Optional: concentration (ng/µL)", min_value=0.0, value=0.0)
    try:
        mw = mb.mw_from_length(int(length), mol, False, mw_method)
        if unit in mb.MASS_UNITS:
            grams = qty * mb.MASS_UNITS[unit]
            moles = grams / mw
        else:
            moles = qty * mb.MOLE_UNITS[unit]
            grams = moles * mw
        st.write(f"MW = **{sci(mw)} g/mol** · copies = **{sci(moles * mb.AVOGADRO)}**")
        cA, cB = st.columns(2)
        cA.dataframe(pd.DataFrame({"Mass": [f"{sci(grams / f)} {u}" for u, f in mb.MASS_UNITS.items()]}),
                     hide_index=True)
        cB.dataframe(pd.DataFrame({"Amount": [f"{sci(moles / f)} {u}" for u, f in mb.MOLE_UNITS.items()]}),
                     hide_index=True)
        if conc > 0:
            nm = mb.ng_per_ul_to_nM(conc, mw)
            st.success(f"{conc:g} ng/µL = {sci(nm)} nM = {sci(nm / 1000)} µM "
                       f"= {sci(mb.copies_from_ng(conc, mw))} copies/µL")
    except ValueError as e:
        st.error(str(e))

# ============================================================================
# 4. SEQUENCE PROPERTIES
# ============================================================================

with tabs[3]:
    st.markdown("### Sequence Properties & Exact Molecular Weight")
    raw = st.text_area("Paste a sequence (plain or FASTA; numbers and spaces are ignored)",
                       height=150, key="sp_seq")
    seq = mb.clean_sequence(raw)
    if seq:
        molecule, problems = mb.validate(seq, allow_ambiguous=False)
        if not show_problems(problems):
            c1, c2 = st.columns(2)
            if molecule == "RNA":
                stype = "RNA"
                c1.info("U detected → treated as RNA.")
                rna_end = c2.radio("5′ end", ["triphosphate (in-vitro transcript)",
                                             "OH (synthetic oligo)"])
                mw_exact = mb.mw_rna_exact(seq, "triphosphate" if rna_end.startswith("tri") else "OH")
            else:
                stype = c1.radio("Treat as", ["dsDNA", "ssDNA"], horizontal=True)
                if stype == "dsDNA":
                    topo = c2.radio("Form", ["Linear, 5′-OH (e.g. PCR product)",
                                             "Linear, 5′-phosphate (e.g. restriction fragment)",
                                             "Circular"])
                    mw_exact = mb.mw_dsdna_exact(seq, five_prime_phosphate="phosphate" in topo,
                                                 circular=topo == "Circular")
                else:
                    p5 = c2.checkbox("5′-phosphate")
                    mw_exact = mb.mw_ssdna_exact(seq, p5)
            L = len(seq)
            unit = LEN_UNIT[stype]
            mw_avg = mb.mw_from_length(L, stype, False, mw_method)
            comp = {b: seq.count(b) for b in ("ACGU" if molecule == "RNA" else "ACGT")}
            m1, m2, m3, m4 = st.columns(4)
            m1.metric("Length", f"{L:,} {unit}")
            m2.metric("GC content", f"{mb.gc_percent(seq):.1f} %")
            m3.metric("Exact MW (from sequence)", f"{sci(mw_exact, 5)} g/mol")
            m4.metric("Length-based MW", f"{sci(mw_avg, 5)} g/mol")
            st.write("Composition: " + " · ".join(f"{b} {n}" for b, n in comp.items()))
            st.write(f"1 µg = **{sci(mb.ng_to_pmol(1000, mw_exact))} pmol**  ·  "
                     f"1 pmol = **{sci(mb.pmol_to_ng(1, mw_exact))} ng**")
            mass = st.number_input("Optional: mass (ng)", min_value=0.0, value=0.0, key="sp_mass")
            if mass > 0:
                st.success(f"{mass:g} ng = {sci(mb.ng_to_pmol(mass, mw_exact))} pmol "
                           f"= {sci(mb.copies_from_ng(mass, mw_exact))} molecules")
            with st.expander("How the exact MW is calculated"):
                st.markdown("""
Single strand: Σ residue masses − 61.96 (5′-OH, 3′-OH), residues
dA 313.21, dC 289.18, dG 329.21, dT 304.20; + 79.98 for a 5′-phosphate.
dsDNA = top strand + reverse-complement strand. Circular DNA has no free ends.
RNA: residues A 329.21, C 305.18, G 345.21, U 306.17; + 159.0 for a 5′-triphosphate
(in-vitro transcript) or − 61.96 for a 5′-OH oligo.
""")

# ============================================================================
# 5. SEQUENCE TOOLS
# ============================================================================

with tabs[4]:
    st.markdown("### Complement / Reverse Complement")
    raw = st.text_area("Sequence (IUPAC codes allowed)", key="rc_seq")
    seq = mb.clean_sequence(raw)
    if seq:
        molecule, problems = mb.validate(seq, allow_ambiguous=True)
        if not show_problems(problems):
            rna = molecule == "RNA"
            st.markdown("**Reverse complement** (5′→3′ — use this for primers/oligos)")
            st.code(wrap(mb.reverse_complement(seq, rna)))
            st.markdown("**Complement** (written 3′→5′, aligned under the input)")
            st.code(wrap(mb.complement(seq, rna)))

    st.markdown("---")
    st.markdown("### Pairwise Alignment")
    amode = st.radio("Mode", ["Primer binding on a template (both strands)",
                              "Compare two sequences (global)"], horizontal=True)
    s1 = mb.clean_sequence(st.text_area("Primer (5′→3′)" if amode.startswith("Primer")
                                        else "Sequence 1", key="al1"))
    s2 = mb.clean_sequence(st.text_area("Template" if amode.startswith("Primer")
                                        else "Sequence 2", key="al2"))
    if amode.startswith("Compare"):
        vs_rc = st.checkbox("Align against the reverse complement of sequence 2 "
                            "(tests antiparallel complementarity)")
    st.caption("Scoring: match +5, mismatch −4, gap open −10, gap extend −0.5 "
               "(EMBOSS needle DNAfull defaults). Affine-gap Needleman–Wunsch/Gotoh.")
    if st.button("Run alignment"):
        probs = []
        for s, name in ((s1, "First sequence"), (s2, "Second sequence")):
            _, p = mb.validate(s)
            probs += [f"{name}: {x}" for x in p]
        if not show_problems(probs):
            try:
                if amode.startswith("Primer"):
                    sites = mb.primer_binding(s1, s2)
                    best = sites[0]
                    ident = 100 * best.aln.matches / len(s1)
                    st.write(f"Best site: **{'forward (matches top strand, anneals to bottom)' if best.strand == '+' else 'reverse (matches bottom strand, anneals to top)'}**, "
                             f"template positions **{best.start}–{best.end}**")
                    m1, m2, m3 = st.columns(3)
                    m1.metric("Primer identity to site", f"{ident:.1f} %")
                    m2.metric("Mismatches + gaps", best.mismatches_total)
                    m3.metric("Errors in 3′-terminal 5 nt", best.mismatches_3p5)
                    if best.terminal_mismatch:
                        st.error("The 3′-terminal base is mismatched — extension will be "
                                 "strongly inhibited (exploited in allele-specific PCR).")
                    elif best.mismatches_3p5:
                        st.warning("Mismatch(es) within the last 5 bases at the 3′ end "
                                   "can reduce or abolish extension.")
                    elif best.mismatches_total:
                        st.info("Mismatches only toward the 5′ end — usually tolerated, "
                                "but they lower the effective Tm.")
                    else:
                        st.success("Perfect match.")
                    st.code(f"Primer   5' {best.aln.a} 3'\n            {best.aln.match_line()}\n"
                            f"Template 5' {best.aln.b} 3'  ({'top' if best.strand == '+' else 'bottom'} strand)")
                    other = sites[1]
                    st.caption(f"Best site on the other strand: score {other.aln.score:.1f} "
                               f"(vs {best.aln.score:.1f}), {other.mismatches_total} mismatches/gaps.")
                else:
                    b = mb.reverse_complement(s2) if vs_rc else s2
                    aln = mb.align(s1, b)
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("Score", f"{aln.score:.1f}")
                    m2.metric("Identity (of alignment length)",
                              f"{100 * aln.matches / aln.length:.1f} %")
                    m3.metric("Identity (of shorter sequence)",
                              f"{100 * aln.matches / min(len(s1), len(s2)):.1f} %")
                    m4.metric("Gaps / mismatches", f"{aln.gaps} / {aln.mismatches}")
                    if vs_rc:
                        st.caption("Sequence 2 was reverse-complemented, so identity here = "
                                   "antiparallel complementarity between the two inputs.")
                    blocks = []
                    for i in range(0, aln.length, 60):
                        blocks.append(f"{aln.a[i:i + 60]}\n{aln.match_line()[i:i + 60]}\n{aln.b[i:i + 60]}")
                    st.code("\n\n".join(blocks))
            except ValueError as e:
                st.error(f"{e} For long sequences use NCBI BLAST or EMBOSS needle.")

# ============================================================================
# 6. PRIMERS
# ============================================================================


def primer_table(rep: mb.PrimerReport, label: str):
    st.markdown(f"**{label}** ({rep.length} nt annealing"
                f"{f' + {len(rep.tail)} nt 5′ tail' if rep.tail else ''})")
    shown = (rep.tail.lower() + rep.anneal) if rep.tail else rep.anneal
    st.code(f"5'-{shown}-3'")
    data = {
        "Tm (nearest-neighbour)": f"{rep.tm:.1f} °C",
        "Tm full length incl. tail": f"{rep.tm_full:.1f} °C" if rep.tm_full else "–",
        "GC content": f"{rep.gc:.0f} %",
        "G/C in last 5 bases": f"{rep.gc_last5}/5" + (" (ends in G/C)" if rep.terminal_gc else ""),
        "Longest mononucleotide run": rep.homopolymer,
        "Longest dinucleotide repeat": rep.dinuc,
        "Hairpin ΔG (kcal/mol)": f"{rep.hairpin.dg:.1f}",
        "Self-dimer ΔG / 3′-end (kcal/mol)": f"{rep.self_dimer.dg:.1f} / {rep.self_dimer_3p.dg:.1f}",
        "Wallace Tm (for comparison only)": f"{mb.tm_wallace(rep.anneal):.0f} °C",
    }
    st.dataframe(pd.DataFrame({"Property": list(data), "Value": [str(v) for v in data.values()]}),
                 hide_index=True)
    if rep.warnings:
        for w in rep.warnings:
            st.warning(w)
    else:
        st.success("No design rule violations.")


with tabs[5]:
    st.caption(f"Tm: SantaLucia 1998 nearest-neighbour with salt and Mg²⁺ correction, "
               f"conditions from the sidebar. ΔG at 37 °C via {mb.thermo_backend()}. "
               f"Lowercase = 5′ tail (not counted in annealing Tm).")
    ptool = st.radio("Tool", ["Design primers from a template", "Analyse existing primers"],
                     horizontal=True)

    if ptool.startswith("Analyse"):
        p1 = mb.clean_sequence(st.text_area("Forward primer (5′→3′)", key="pa1"))
        p2 = mb.clean_sequence(st.text_area("Reverse primer (5′→3′, optional)", key="pa2"))
        if st.button("Analyse"):
            probs = []
            for p, nme in ((p1, "Forward"), (p2, "Reverse")):
                if p and (set(p) - set("ACGT")):
                    probs.append(f"{nme} primer: only A, C, G, T allowed for Tm/ΔG.")
            if not p1:
                probs.append("Enter at least one primer.")
            if not show_problems(probs):
                reps = [mb.analyse_primer(p1, cond)]
                cols = st.columns(2 if p2 else 1)
                with cols[0]:
                    primer_table(reps[0], "Forward primer")
                if p2:
                    reps.append(mb.analyse_primer(p2, cond))
                    with cols[1]:
                        primer_table(reps[1], "Reverse primer")
                    het_any, het_3p = mb.dimer(p1, p2, cond)
                    st.markdown("**Pair**")
                    st.write(f"Tm difference: **{abs(reps[0].tm - reps[1].tm):.1f} °C** · "
                             f"hetero-dimer ΔG **{het_any.dg:.1f}** kcal/mol · "
                             f"3′-end dimer ΔG **{het_3p.dg:.1f}** kcal/mol")
                    if het_3p.dg < mb.DIMER_3P_LIMIT or het_any.dg < mb.DIMER_STRONG:
                        st.warning("Primer-dimer risk.")
                st.info("Annealing temperature: for Taq-type enzymes start ~3–5 °C below the "
                        "lower primer Tm (or run a gradient). High-fidelity enzymes "
                        "(Q5, Phusion) use buffer-specific rules — check the supplier's calculator.")

    else:
        tpl = mb.clean_sequence(st.text_area("Template (5′→3′ top strand, plain or FASTA)",
                                             height=150, key="pd_tpl"))
        c1, c2, c3 = st.columns(3)
        with c1:
            dmode = st.radio("Primer placement",
                             ["Anchored at the ends (amplify exactly this sequence, e.g. ORF)",
                              "Anywhere within flanking regions"])
            region = st.number_input("Flanking region (bp)", min_value=20, value=100, step=10,
                                     disabled=dmode.startswith("Anchored"))
        with c2:
            target = st.number_input("Target Tm (°C)", 45.0, 75.0, 60.0, 0.5)
            lmin, lmax = st.slider("Annealing length (nt)", 12, 40, (18, 30))
            maxdiff = st.number_input("Max Tm difference (°C)", 0.5, 10.0, 3.0, 0.5)
        with c3:
            tail_kind = st.selectbox("5′ additions",
                                     ["None", "Restriction sites", "Custom tails (e.g. Gibson overlaps)"])
            tail_f = tail_r = ""
            enz_f = enz_r = None
            if tail_kind == "Restriction sites":
                enz_f = st.selectbox("Forward-primer enzyme", list(mb.RESTRICTION_SITES), index=2)
                enz_r = st.selectbox("Reverse-primer enzyme", list(mb.RESTRICTION_SITES), index=15)
                leader = mb.clean_sequence(st.text_input(
                    "Leader bases 5′ of the site", "TAAGCA",
                    help="Most enzymes cut poorly at the very end of DNA; 4–6 extra bases "
                         "are usually enough (see the NEB 'cleavage close to the end' chart)."))
                tail_f = leader + mb.RESTRICTION_SITES[enz_f]
                tail_r = leader + mb.RESTRICTION_SITES[enz_r]
            elif tail_kind.startswith("Custom"):
                tail_f = mb.clean_sequence(st.text_input("Forward 5′ tail", ""))
                tail_r = mb.clean_sequence(st.text_input("Reverse 5′ tail", ""))

        if st.button("Design primers"):
            probs = []
            if not tpl:
                probs.append("Paste a template.")
            elif set(tpl) - set("ACGT"):
                probs.append("Template must contain only A, C, G, T for primer design.")
            if set(tail_f + tail_r) - set("ACGT"):
                probs.append("Tails must contain only A, C, G, T.")
            if not show_problems(probs):
                try:
                    with st.spinner("Scoring candidate primers…"):
                        pairs = mb.design_primers(
                            tpl, cond, target_tm=target, min_len=lmin, max_len=lmax,
                            mode="anchored" if dmode.startswith("Anchored") else "flanking",
                            region=int(region), max_tm_diff=maxdiff,
                            tail_f=tail_f, tail_r=tail_r)
                    if not pairs:
                        st.error("No valid primer pair found — widen the length range or region.")
                    for enz in {e for e in (enz_f, enz_r) if e}:
                        hits = mb.find_sites(tpl, mb.RESTRICTION_SITES[enz])
                        if hits:
                            st.error(f"{enz} also cuts inside the template at position(s) "
                                     f"{', '.join(map(str, hits))} — choose another enzyme.")
                    for k, p in enumerate(pairs, 1):
                        with st.expander(
                                f"Pair {k}: product {p.product_bp} bp · Tm {p.fwd.tm:.1f} / "
                                f"{p.rev.tm:.1f} °C · penalty {p.penalty:.1f}", expanded=k == 1):
                            a, b = st.columns(2)
                            with a:
                                primer_table(p.fwd, f"Forward (template {p.fwd_start + 1}–"
                                                    f"{p.fwd_start + p.fwd.length})")
                            with b:
                                primer_table(p.rev, f"Reverse (template {p.rev_end - p.rev.length + 1}–"
                                                    f"{p.rev_end})")
                            st.write(f"Hetero-dimer ΔG {p.hetero.dg:.1f} kcal/mol · "
                                     f"3′-end {p.hetero_3p.dg:.1f} kcal/mol · "
                                     f"Tm difference {p.tm_diff:.1f} °C")
                            for w in p.warnings:
                                st.warning(w)
                    if pairs:
                        st.info("Check specificity against the genome/plasmid background "
                                "(e.g. NCBI Primer-BLAST) before ordering. In anchored mode the "
                                "primer start is fixed, so GC-rich or AT-rich ends may force "
                                "rule violations — consider flanking mode or a longer/shorter range.")
                except ValueError as e:
                    st.error(str(e))

st.markdown("---")
st.caption("🧬 RNMolBioTools • Molecular Biology Utilities · V9 — calculations are planning "
           "aids; verify critical values experimentally.")
