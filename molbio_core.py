"""
molbio_core.py — scientific core of RNMolBioTools (DNA Cloning Calculator)

All calculations live here so they can be tested independently of the GUI.
No third-party dependency is required. If `primer3-py` is installed, its
thermodynamic hairpin/dimer calculations are used; otherwise a built-in
nearest-neighbour approximation is used.

Conventions (stated once, used everywhere)
------------------------------------------
Molecular weights (anhydrous free acid, as used by NEB / QIAGEN calculators):
    dsDNA (linear)   MW = L(bp) * 617.96 + 36.04   g/mol
    dsDNA (circular) MW = L(bp) * 617.96           (no free ends)
    ssDNA            MW = L(nt) * 308.97 + 18.02
    ssRNA            MW = L(nt) * 321.47 + 18.02
Legacy averages (sodium salt, "rule of thumb"): 650 g/mol/bp, 330 g/mol/nt
(ssDNA), 340 g/mol/nt (RNA) — available as an option only.

Melting temperature: SantaLucia (1998) unified nearest-neighbour parameters,
SantaLucia salt correction of entropy, Mg2+ converted to Na+-equivalent with
the von Ahsen et al. (2001) relation  [Na+]eq = [Na+] + 120*sqrt([Mg2+]-[dNTP])
(the same approach Primer3 uses by default). Primer is assumed to be in large
excess over template, so Tm = dH / (dS + R ln C_primer).
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

AVOGADRO = 6.02214076e23
R_CAL = 1.9872  # cal / (K mol)

# --------------------------------------------------------------------------- #
# Sequence handling
# --------------------------------------------------------------------------- #

IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")

_COMP_DNA = {
    "A": "T", "T": "A", "G": "C", "C": "G", "U": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D", "N": "N", "-": "-",
}


def clean_sequence(raw: str) -> str:
    """Strip FASTA header lines, GenBank-style numbers, whitespace; uppercase."""
    if not raw:
        return ""
    lines = [ln for ln in raw.splitlines() if not ln.strip().startswith(">")]
    seq = "".join(lines)
    seq = re.sub(r"[\s\d/]", "", seq)
    return seq.upper()


def validate(seq: str, allow_ambiguous: bool = True) -> tuple[str, list[str]]:
    """Return (molecule, problems). molecule is 'DNA', 'RNA' or '?'."""
    problems = []
    if not seq:
        return "?", ["Sequence is empty."]
    has_t, has_u = "T" in seq, "U" in seq
    if has_t and has_u:
        problems.append("Sequence contains both T and U.")
    molecule = "RNA" if has_u and not has_t else "DNA"
    alphabet = IUPAC_RNA if molecule == "RNA" else IUPAC_DNA
    bad = sorted(set(seq) - alphabet)
    if bad:
        problems.append(f"Invalid characters: {', '.join(bad)}")
    if not allow_ambiguous:
        strict = set("ACGU") if molecule == "RNA" else set("ACGT")
        amb = sorted(set(seq) - strict - set(bad))
        if amb:
            problems.append(
                f"Ambiguous IUPAC codes ({', '.join(amb)}) cannot be used for this calculation."
            )
    return molecule, problems


def complement(seq: str, rna: bool = False) -> str:
    out = "".join(_COMP_DNA.get(b, "?") for b in seq)
    return out.replace("T", "U") if rna else out


def reverse_complement(seq: str, rna: bool = False) -> str:
    return complement(seq, rna)[::-1]


def gc_percent(seq: str) -> float:
    if not seq:
        return 0.0
    return 100.0 * (seq.count("G") + seq.count("C") + seq.count("S")) / len(seq)


# --------------------------------------------------------------------------- #
# Molecular weight
# --------------------------------------------------------------------------- #

MW_METHODS = ("Exact (NEB/QIAGEN formula)", "Legacy average (650/330/340)")


def mw_from_length(length: int, molecule: str = "dsDNA", circular: bool = False,
                   method: str = MW_METHODS[0]) -> float:
    """Molecular weight (g/mol) from length only. molecule: dsDNA | ssDNA | RNA."""
    if length <= 0:
        raise ValueError("Length must be > 0.")
    legacy = method == MW_METHODS[1]
    if molecule == "dsDNA":
        if legacy:
            return length * 650.0
        return length * 617.96 + (0.0 if circular else 36.04)
    if molecule == "ssDNA":
        return length * (330.0 if legacy else 308.97) + (0.0 if legacy else 18.02)
    if molecule == "RNA":
        return length * (340.0 if legacy else 321.47) + (0.0 if legacy else 18.02)
    raise ValueError(f"Unknown molecule type: {molecule}")


# Anhydrous nucleotide residue masses (g/mol)
_DNA_RES = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.20}
_RNA_RES = {"A": 329.21, "C": 305.18, "G": 345.21, "U": 306.17}
PHOSPHATE = 79.98  # HPO3 added for a 5'-monophosphate


def mw_ssdna_exact(seq: str, five_prime_phosphate: bool = False) -> float:
    """Exact MW of a single DNA strand (5'-OH, 3'-OH unless 5'-phosphate)."""
    mw = sum(_DNA_RES[b] for b in seq) - 61.96
    return mw + (PHOSPHATE if five_prime_phosphate else 0.0)


def mw_rna_exact(seq: str, five_prime: str = "triphosphate") -> float:
    """Exact RNA MW. five_prime: 'triphosphate' (in-vitro transcript) or 'OH' (synthetic)."""
    s = sum(_RNA_RES[b] for b in seq)
    return s + 159.0 if five_prime == "triphosphate" else s - 61.96


def mw_dsdna_exact(seq: str, five_prime_phosphate: bool = False,
                   circular: bool = False) -> float:
    """Exact MW of a dsDNA duplex = top strand + bottom strand."""
    top = mw_ssdna_exact(seq)
    bottom = mw_ssdna_exact(reverse_complement(seq))
    mw = top + bottom
    if circular:
        # closing each strand removes one H2O (18.02) relative to the 5'-OH/3'-OH form
        # and adds one phosphodiester (HPO3, 79.98)
        mw += 2 * (PHOSPHATE - 18.02)
    elif five_prime_phosphate:
        mw += 2 * PHOSPHATE
    return mw


# --------------------------------------------------------------------------- #
# Mass / moles / copies
# --------------------------------------------------------------------------- #

MASS_UNITS = {"g": 1.0, "mg": 1e-3, "µg": 1e-6, "ng": 1e-9, "pg": 1e-12, "fg": 1e-15}
MOLE_UNITS = {"mol": 1.0, "mmol": 1e-3, "µmol": 1e-6, "nmol": 1e-9,
              "pmol": 1e-12, "fmol": 1e-15, "amol": 1e-18}


def ng_to_pmol(ng: float, mw: float) -> float:
    return ng * 1e3 / mw


def pmol_to_ng(pmol: float, mw: float) -> float:
    return pmol * mw / 1e3


def copies_from_ng(ng: float, mw: float) -> float:
    return ng * 1e-9 / mw * AVOGADRO


def ng_per_ul_to_nM(conc_ng_ul: float, mw: float) -> float:
    return conc_ng_ul * 1e6 / mw


def haploid_genome_mass_pg(genome_bp: float) -> float:
    return genome_bp * 617.96 / AVOGADRO * 1e12


# --------------------------------------------------------------------------- #
# Ligation / assembly
# --------------------------------------------------------------------------- #

@dataclass
class Component:
    name: str
    bp: int
    ratio: float          # molar parts relative to vector = 1
    mw: float = 0.0
    ng: float = 0.0
    pmol: float = 0.0
    stock_ng_ul: float | None = None
    ul: float | None = None


@dataclass
class LigationResult:
    components: list
    total_ng: float
    total_pmol: float
    dna_volume_ul: float | None
    warnings: list = field(default_factory=list)


def ligation_mix(vector_bp: int, inserts: list[tuple[int, float]], mode: str, amount: float,
                 vector_circular_before_cut: bool = False, stocks: list | None = None,
                 reaction_ul: float | None = None, other_ul: float = 0.0) -> LigationResult:
    """
    inserts: list of (insert_bp, molar ratio insert:vector), e.g. [(1000, 3)].
    mode: 'vector_ng'  -> amount = vector mass (ng) [recommended]
          'total_ng'   -> amount = total DNA mass (ng)
          'total_pmol' -> amount = total pmol of all fragments (assembly kits)
    stocks: optional list of stock concentrations (ng/µL), vector first.
    """
    if vector_bp <= 0 or any(bp <= 0 for bp, _ in inserts):
        raise ValueError("All fragment sizes must be > 0 bp.")
    if any(r <= 0 for _, r in inserts):
        raise ValueError("Molar ratios must be > 0.")
    if amount <= 0:
        raise ValueError("Amount must be > 0.")

    vec = Component("Vector", vector_bp, 1.0, mw_from_length(vector_bp, "dsDNA"))
    ins = [Component(f"Insert {k + 1}" if len(inserts) > 1 else "Insert", bp, r,
                     mw_from_length(bp, "dsDNA")) for k, (bp, r) in enumerate(inserts)]

    if mode == "vector_ng":
        vec.ng = amount
        vec.pmol = ng_to_pmol(vec.ng, vec.mw)
    elif mode == "total_ng":
        # total = v_ng * (1 + sum r_i * MW_i / MW_v)
        vec.ng = amount / (1 + sum(c.ratio * c.mw / vec.mw for c in ins))
        vec.pmol = ng_to_pmol(vec.ng, vec.mw)
    elif mode == "total_pmol":
        vec.pmol = amount / (1 + sum(c.ratio for c in ins))
        vec.ng = pmol_to_ng(vec.pmol, vec.mw)
    else:
        raise ValueError(mode)

    for c in ins:
        c.pmol = c.ratio * vec.pmol
        c.ng = pmol_to_ng(c.pmol, c.mw)

    comps = [vec] + ins
    warnings = []
    dna_ul = None
    if stocks and len(stocks) == len(comps) and all(s and s > 0 for s in stocks):
        for c, s in zip(comps, stocks):
            c.stock_ng_ul = s
            c.ul = c.ng / s
            if c.ul < 0.5:
                warnings.append(f"{c.name}: {c.ul:.2f} µL is below what can be pipetted "
                                f"accurately — dilute the stock first.")
        dna_ul = sum(c.ul for c in comps)
        if reaction_ul:
            free = reaction_ul - other_ul
            if dna_ul > free:
                warnings.append(f"DNA volume ({dna_ul:.2f} µL) exceeds the space left in the "
                                f"reaction ({free:.2f} µL). Use more concentrated stocks or "
                                f"a larger reaction volume.")
    total_ng = sum(c.ng for c in comps)
    total_pmol = sum(c.pmol for c in comps)
    if reaction_ul and total_ng / reaction_ul > 10:
        warnings.append(f"Total DNA is {total_ng / reaction_ul:.1f} ng/µL; very high DNA "
                        f"concentrations favour concatemers over circularisation.")
    return LigationResult(comps, total_ng, total_pmol, dna_ul, warnings)


# --------------------------------------------------------------------------- #
# Serial dilution (qPCR standards)
# --------------------------------------------------------------------------- #

def serial_dilution(stock_conc: float, top_conc: float, fold: float, n_points: int,
                    final_ul: float) -> tuple[list[dict], list[str]]:
    """
    Each tube keeps `final_ul` after transferring to the next tube.
    Returns rows with transfer volume, diluent volume and concentration.
    """
    warnings = []
    if top_conc > stock_conc:
        raise ValueError("Top standard cannot be more concentrated than the stock.")
    if fold <= 1 or n_points < 1 or final_ul <= 0:
        raise ValueError("Fold must be > 1, points ≥ 1, volume > 0.")
    totals = [0.0] * n_points
    totals[-1] = final_ul
    for i in range(n_points - 2, -1, -1):
        totals[i] = final_ul + totals[i + 1] / fold
    rows = []
    for i in range(n_points):
        conc = top_conc / fold ** i
        if i == 0:
            transfer = totals[0] * top_conc / stock_conc
            source = "stock"
        else:
            transfer = totals[i] / fold
            source = f"tube {i}"
        rows.append({"Tube": i + 1, "Concentration": conc, "Transfer (µL)": transfer,
                     "From": source, "Diluent (µL)": totals[i] - transfer,
                     "Total made (µL)": totals[i]})
    if rows[0]["Transfer (µL)"] < 1:
        warnings.append(f"Tube 1 needs {rows[0]['Transfer (µL)']:.3g} µL of stock — make an "
                        f"intermediate pre-dilution of the stock first (e.g. 1:100).")
    return rows, warnings


# --------------------------------------------------------------------------- #
# Nearest-neighbour thermodynamics (SantaLucia 1998 unified)
# --------------------------------------------------------------------------- #

# dH (kcal/mol), dS (cal/K/mol); key = 5'->3' top-strand dinucleotide
_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_INIT = {"GC": (0.1, -2.8), "AT": (2.3, 4.1)}   # per terminal base pair
_SYM = (0.0, -1.4)

# Hairpin loop free energies at 37 °C (SantaLucia & Hicks 2004), kcal/mol
_HP_LOOP = {3: 3.5, 4: 3.5, 5: 3.3, 6: 4.0, 7: 4.2, 8: 4.3, 9: 4.5, 10: 4.6,
            12: 5.0, 14: 5.1, 16: 5.3, 18: 5.5, 20: 5.7, 25: 6.1, 30: 6.3}


@dataclass
class PCRConditions:
    na_mM: float = 50.0       # monovalent cations
    mg_mM: float = 1.5        # free + bound Mg2+
    dntp_mM: float = 0.8      # total dNTP (0.2 mM each)
    primer_nM: float = 250.0  # each primer

    def na_equivalent_M(self) -> float:
        free_mg = max(self.mg_mM - self.dntp_mM, 0.0)
        return (self.na_mM + 120.0 * math.sqrt(free_mg)) / 1000.0


def _terminal(b: str) -> str:
    return "GC" if b in "GC" else "AT"


def _duplex_dh_ds(seq: str) -> tuple[float, float]:
    """dH, dS of a perfectly matched duplex of `seq` with its complement (1 M Na+)."""
    dh = ds = 0.0
    for i in range(len(seq) - 1):
        h, s = _NN[seq[i:i + 2]]
        dh += h
        ds += s
    for b in (seq[0], seq[-1]):
        h, s = _INIT[_terminal(b)]
        dh += h
        ds += s
    if seq == reverse_complement(seq):
        dh += _SYM[0]
        ds += _SYM[1]
    return dh, ds


def _salt_corrected_ds(ds: float, n_bp: int, cond: PCRConditions) -> float:
    return ds + 0.368 * (n_bp - 1) * math.log(cond.na_equivalent_M())


def tm_nn(seq: str, cond: PCRConditions | None = None) -> float:
    """Nearest-neighbour Tm (°C) of a primer against its perfect complement."""
    cond = cond or PCRConditions()
    if len(seq) < 2 or set(seq) - set("ACGT"):
        raise ValueError("Tm needs an unambiguous A/C/G/T sequence of ≥ 2 nt.")
    dh, ds = _duplex_dh_ds(seq)
    ds = _salt_corrected_ds(ds, len(seq), cond)
    c = cond.primer_nM * 1e-9   # primer in large excess over template -> ln(C)
    return dh * 1000.0 / (ds + R_CAL * math.log(c)) - 273.15


def tm_wallace(seq: str) -> float:
    """Wallace rule — shown for comparison only (valid ~14–20 nt, high salt)."""
    return 2 * (seq.count("A") + seq.count("T")) + 4 * (seq.count("G") + seq.count("C"))


def _stack_dg(seq_top: str, cond: PCRConditions, t_c: float = 37.0) -> float:
    """dG of contiguous WC stacks for the top-strand segment (no initiation)."""
    t = t_c + 273.15
    dh = ds = 0.0
    for i in range(len(seq_top) - 1):
        h, s = _NN[seq_top[i:i + 2]]
        dh += h
        ds += s
    ds = _salt_corrected_ds(ds, len(seq_top), cond) if len(seq_top) > 1 else ds
    return dh - t * ds / 1000.0


def _init_dg(b1: str, b2: str, t_c: float = 37.0) -> float:
    t = t_c + 273.15
    g = 0.0
    for b in (b1, b2):
        h, s = _INIT[_terminal(b)]
        g += h - t * s / 1000.0
    return g


def _loop_dg(n: int, t_c: float = 37.0) -> float:
    if n in _HP_LOOP:
        return _HP_LOOP[n]
    keys = sorted(_HP_LOOP)
    if n > keys[-1]:
        ref = keys[-1]
    else:
        ref = max(k for k in keys if k < n)
    return _HP_LOOP[ref] + 2.44 * R_CAL * (t_c + 273.15) * math.log(n / ref) / 1000.0


@dataclass
class Structure:
    dg: float                 # kcal/mol at 37 °C (more negative = more stable)
    three_prime: bool         # does the structure involve a primer 3' end?
    description: str


def dimer_dg(a: str, b: str, cond: PCRConditions | None = None,
             min_bp: int = 3) -> tuple[Structure, Structure]:
    """
    Built-in ungapped duplex search between primers a and b (both 5'->3').
    Returns (most stable duplex anywhere, most stable duplex involving a 3' end).
    """
    cond = cond or PCRConditions()
    brev = b[::-1]          # b written 3'->5' so it lies antiparallel under a
    la, lb = len(a), len(brev)
    best_any = Structure(0.0, False, "none")
    best_3p = Structure(0.0, True, "none")
    for k in range(-(lb - 1), la):
        run_start = None
        for i in range(max(0, k), min(la, k + lb) + 1):
            paired = i < min(la, k + lb) and _COMP_DNA.get(a[i]) == brev[i - k]
            if paired and run_start is None:
                run_start = i
            if (not paired or i == min(la, k + lb)) and run_start is not None:
                end = i  # exclusive
                n = end - run_start
                if n >= min_bp:
                    seg = a[run_start:end]
                    dg = _stack_dg(seg, cond) + _init_dg(seg[0], seg[-1])
                    a3 = end == la                       # a's 3' terminal base paired
                    b3 = (run_start - k) == 0            # b's 3' terminal base paired
                    desc = f"{n} bp: 5'-{seg}-3'"
                    if dg < best_any.dg:
                        best_any = Structure(dg, a3 or b3, desc)
                    if (a3 or b3) and dg < best_3p.dg:
                        best_3p = Structure(dg, True, desc)
                run_start = None
    return best_any, best_3p


def hairpin_dg(s: str, cond: PCRConditions | None = None, min_loop: int = 3,
               min_stem: int = 3) -> Structure:
    """Built-in hairpin search: perfect stems ≥ min_stem bp, loops ≥ min_loop nt."""
    cond = cond or PCRConditions()
    best = Structure(0.0, False, "none")
    n = len(s)
    for i in range(n):
        for j in range(n - 1, i + min_loop + 2 * min_stem - 2, -1):
            k = 0
            while (i + k < j - k and _COMP_DNA.get(s[i + k]) == s[j - k]
                   and (j - k) - (i + k) - 1 >= min_loop):
                k += 1
                if k >= min_stem:
                    loop = (j - k + 1) - (i + k)
                    dg = _stack_dg(s[i:i + k], cond) + _loop_dg(loop)
                    if dg < best.dg:
                        best = Structure(dg, j == n - 1,
                                         f"{k} bp stem, {loop} nt loop")
            # continue scanning
    return best


# Optional primer3-py backend ------------------------------------------------ #
try:  # pragma: no cover - depends on environment
    import primer3 as _p3
    _HAS_PRIMER3 = True
except Exception:  # noqa: BLE001
    _p3 = None
    _HAS_PRIMER3 = False


def thermo_backend() -> str:
    return "primer3-py" if _HAS_PRIMER3 else "built-in nearest-neighbour approximation"


def _p3_kwargs(cond: PCRConditions) -> dict:
    return dict(mv_conc=cond.na_mM, dv_conc=cond.mg_mM, dntp_conc=cond.dntp_mM,
                dna_conc=cond.primer_nM)


def hairpin(seq: str, cond: PCRConditions) -> Structure:
    if _HAS_PRIMER3:
        try:
            r = _p3.calc_hairpin(seq, **_p3_kwargs(cond))
            return Structure(r.dg / 1000.0 if r.structure_found else 0.0, False,
                             "primer3 thermodynamic model")
        except Exception:  # noqa: BLE001
            pass
    return hairpin_dg(seq, cond)


def dimer(a: str, b: str, cond: PCRConditions) -> tuple[Structure, Structure]:
    builtin_any, builtin_3p = dimer_dg(a, b, cond)
    if _HAS_PRIMER3:
        try:
            r = (_p3.calc_homodimer(a, **_p3_kwargs(cond)) if a == b
                 else _p3.calc_heterodimer(a, b, **_p3_kwargs(cond)))
            any_dg = r.dg / 1000.0 if r.structure_found else 0.0
            e1 = _p3.calc_end_stability(a, b, **_p3_kwargs(cond)).dg / 1000.0
            e2 = _p3.calc_end_stability(b, a, **_p3_kwargs(cond)).dg / 1000.0
            return (Structure(any_dg, False, "primer3 thermodynamic model"),
                    Structure(min(e1, e2, 0.0), True, "primer3 3'-end stability"))
        except Exception:  # noqa: BLE001
            pass
    return builtin_any, builtin_3p


# Thresholds (kcal/mol) — widely used primer-design guidance
HAIRPIN_3P_LIMIT = -2.0
HAIRPIN_INT_LIMIT = -3.0
DIMER_3P_LIMIT = -5.0
DIMER_INT_LIMIT = -6.0
DIMER_STRONG = -9.0


# --------------------------------------------------------------------------- #
# Primer metrics & design
# --------------------------------------------------------------------------- #

def max_homopolymer(seq: str) -> int:
    best = run = 1 if seq else 0
    for i in range(1, len(seq)):
        run = run + 1 if seq[i] == seq[i - 1] else 1
        best = max(best, run)
    return best


def max_dinuc_repeat(seq: str) -> int:
    """Longest run of a repeated dinucleotide (e.g. ATATAT = 3)."""
    best = 0
    for start in (0, 1):
        i = start
        while i + 1 < len(seq):
            unit = seq[i:i + 2]
            if unit[0] == unit[1]:
                i += 1
                continue
            n = 1
            while seq[i + 2 * n:i + 2 * n + 2] == unit:
                n += 1
            best = max(best, n)
            i += 2 * n if n > 1 else 1
    return best


@dataclass
class PrimerReport:
    anneal: str
    tail: str
    length: int
    gc: float
    tm: float
    tm_full: float | None
    gc_last5: int
    terminal_gc: bool
    homopolymer: int
    dinuc: int
    hairpin: Structure
    self_dimer: Structure
    self_dimer_3p: Structure
    warnings: list
    penalty: float

    @property
    def full(self) -> str:
        return self.tail + self.anneal


def analyse_primer(anneal: str, cond: PCRConditions, tail: str = "",
                   target_tm: float | None = None) -> PrimerReport:
    full = tail + anneal
    tm = tm_nn(anneal, cond)
    tm_full = tm_nn(full, cond) if tail else None
    gc = gc_percent(anneal)
    last5 = anneal[-5:]
    gc5 = last5.count("G") + last5.count("C")
    term = anneal[-1] in "GC"
    hp = max_homopolymer(anneal)
    dn = max_dinuc_repeat(anneal)
    hpin = hairpin(full, cond)
    sd_any, sd_3p = dimer(full, full, cond)

    w = []
    pen = 0.0
    if not 40 <= gc <= 60:
        w.append(f"GC content {gc:.0f}% is outside 40–60%.")
        pen += 0.2 * min(abs(gc - 40), abs(gc - 60))
    if gc5 == 0:
        w.append("No G/C in the last 5 bases (weak 3′ clamp).")
        pen += 2
    elif gc5 > 3:
        w.append(f"{gc5}/5 G/C in the last 5 bases — too GC-rich 3′ end, risk of mispriming.")
        pen += 2 * (gc5 - 3)
    if hp >= 5:
        w.append(f"Run of {hp} identical bases (keep ≤ 4).")
        pen += 3
    if dn >= 5:
        w.append(f"Dinucleotide repeat ×{dn} (keep ≤ 4).")
        pen += 3
    hp_limit = HAIRPIN_3P_LIMIT if hpin.three_prime else HAIRPIN_INT_LIMIT
    if hpin.dg < hp_limit:
        w.append(f"Hairpin ΔG {hpin.dg:.1f} kcal/mol ({hpin.description})"
                 f"{' involving the 3′ end' if hpin.three_prime else ''}.")
        pen += hp_limit - hpin.dg
    if sd_3p.dg < DIMER_3P_LIMIT:
        w.append(f"Self-dimer at the 3′ end, ΔG {sd_3p.dg:.1f} kcal/mol.")
        pen += DIMER_3P_LIMIT - sd_3p.dg
    if sd_any.dg < DIMER_INT_LIMIT:
        w.append(f"Self-dimer ΔG {sd_any.dg:.1f} kcal/mol.")
        pen += 0.5 * (DIMER_INT_LIMIT - sd_any.dg)
    if target_tm is not None:
        pen += abs(tm - target_tm)
    return PrimerReport(anneal, tail, len(anneal), gc, tm, tm_full, gc5, term, hp, dn,
                        hpin, sd_any, sd_3p, w, pen)


@dataclass
class PairReport:
    fwd: PrimerReport
    rev: PrimerReport
    fwd_start: int          # 0-based on template
    rev_end: int            # exclusive, 0-based on template
    product_bp: int
    tm_diff: float
    hetero: Structure
    hetero_3p: Structure
    warnings: list
    penalty: float


def pair_primers(f: PrimerReport, r: PrimerReport, fs: int, re_: int,
                 cond: PCRConditions, max_tm_diff: float) -> PairReport:
    het_any, het_3p = dimer(f.full, r.full, cond)
    w = []
    d = abs(f.tm - r.tm)
    pen = f.penalty + r.penalty + 1.5 * d
    if d > max_tm_diff:
        w.append(f"Tm difference {d:.1f} °C exceeds {max_tm_diff:.0f} °C.")
        pen += 5 * (d - max_tm_diff)
    if het_3p.dg < DIMER_3P_LIMIT:
        w.append(f"Primer-dimer at a 3′ end, ΔG {het_3p.dg:.1f} kcal/mol.")
        pen += DIMER_3P_LIMIT - het_3p.dg
    if het_any.dg < DIMER_INT_LIMIT:
        w.append(f"Primer-dimer ΔG {het_any.dg:.1f} kcal/mol.")
        pen += 0.5 * (DIMER_INT_LIMIT - het_any.dg)
    product = (re_ - fs) + len(f.tail) + len(r.tail)
    return PairReport(f, r, fs, re_, product, d, het_any, het_3p, w, pen)


def design_primers(template: str, cond: PCRConditions, target_tm: float = 60.0,
                   min_len: int = 18, max_len: int = 30, mode: str = "anchored",
                   region: int = 100, max_tm_diff: float = 3.0,
                   tail_f: str = "", tail_r: str = "", n_results: int = 3) -> list[PairReport]:
    """
    mode 'anchored' : primers start exactly at the template ends (amplify the whole
                      sequence, e.g. an ORF for cloning); only the length varies.
    mode 'flanking' : primers may start anywhere within `region` bp of each end.
    """
    n = len(template)
    if set(template) - set("ACGT"):
        raise ValueError("Template must contain only A, C, G, T for primer design.")
    if n < 2 * min_len:
        raise ValueError("Template is too short for two primers.")
    region = min(region, n // 2)

    def candidates(strand_seq: str, starts) -> list:
        out = []
        for s in starts:
            for L in range(min_len, max_len + 1):
                if s + L > len(strand_seq):
                    break
                seq = strand_seq[s:s + L]
                try:
                    tm = tm_nn(seq, cond)
                except ValueError:
                    continue
                if abs(tm - target_tm) <= 6 or mode == "anchored":
                    out.append((s, seq, tm))
        out.sort(key=lambda t: abs(t[2] - target_tm))
        return out[:60]

    rc = reverse_complement(template)
    starts = [0] if mode == "anchored" else range(0, max(1, region - min_len + 1))
    f_cand = candidates(template, starts)
    r_cand = candidates(rc, starts)

    f_rep = sorted(((s, analyse_primer(q, cond, tail_f, target_tm)) for s, q, _ in f_cand),
                   key=lambda t: t[1].penalty)[:15]
    r_rep = sorted(((s, analyse_primer(q, cond, tail_r, target_tm)) for s, q, _ in r_cand),
                   key=lambda t: t[1].penalty)[:15]
    pairs = []
    for fs, f in f_rep:
        for rs, r in r_rep:
            rev_end = n - rs
            if rev_end - fs < f.length + r.length:
                continue
            pairs.append(pair_primers(f, r, fs, rev_end, cond, max_tm_diff))
    pairs.sort(key=lambda p: p.penalty)
    unique, seen = [], set()
    for p in pairs:
        key = (p.fwd.anneal, p.rev.anneal)
        if key not in seen:
            seen.add(key)
            unique.append(p)
        if len(unique) == n_results:
            break
    return unique


RESTRICTION_SITES = {
    "BamHI": "GGATCC", "BglII": "AGATCT", "EcoRI": "GAATTC", "EcoRV": "GATATC",
    "HindIII": "AAGCTT", "KpnI": "GGTACC", "NcoI": "CCATGG", "NdeI": "CATATG",
    "NheI": "GCTAGC", "NotI": "GCGGCCGC", "PstI": "CTGCAG", "SacI": "GAGCTC",
    "SalI": "GTCGAC", "SpeI": "ACTAGT", "XbaI": "TCTAGA", "XhoI": "CTCGAG",
}


def find_sites(template: str, site: str) -> list[int]:
    """1-based positions of a site on either strand."""
    hits = set()
    for m in re.finditer(f"(?={site})", template):
        hits.add(m.start() + 1)
    rc_site = reverse_complement(site)
    if rc_site != site:
        for m in re.finditer(f"(?={rc_site})", template):
            hits.add(m.start() + 1)
    return sorted(hits)


# --------------------------------------------------------------------------- #
# Alignment (Gotoh affine gaps; global or semi-global)
# --------------------------------------------------------------------------- #

MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 5.0, -4.0, 10.0, 0.5   # EMBOSS needle DNAfull
_NEG = -1e18
MAX_CELLS = 1_000_000   # pure-Python DP: keeps memory < ~150 MB on Streamlit Cloud


@dataclass
class Alignment:
    a: str
    b: str
    score: float
    b_start: int   # 0-based start on b (semi-global)
    b_end: int     # exclusive

    @property
    def length(self) -> int:
        return len(self.a)

    @property
    def matches(self) -> int:
        return sum(1 for x, y in zip(self.a, self.b) if x == y and x != "-")

    @property
    def gaps(self) -> int:
        return sum(1 for x, y in zip(self.a, self.b) if x == "-" or y == "-")

    @property
    def mismatches(self) -> int:
        return self.length - self.matches - self.gaps

    def match_line(self) -> str:
        return "".join("|" if x == y and x != "-" else (" " if "-" in (x, y) else ".")
                       for x, y in zip(self.a, self.b))


def align(a: str, b: str, semiglobal: bool = False, match: float = MATCH,
          mismatch: float = MISMATCH, gap_open: float = GAP_OPEN,
          gap_extend: float = GAP_EXTEND) -> Alignment:
    """
    Gotoh affine-gap alignment; gap of length L costs gap_open + (L-1)*gap_extend.
    semiglobal=True: `a` (e.g. primer) must align end-to-end, leading/trailing parts
    of `b` (template) are free.
    """
    n, m = len(a), len(b)
    if n == 0 or m == 0:
        raise ValueError("Both sequences must be non-empty.")
    if n * m > MAX_CELLS:
        raise ValueError(f"Sequences too long for the in-browser aligner "
                         f"({n}×{m} > {MAX_CELLS:,} cells).")
    M = [[_NEG] * (m + 1) for _ in range(n + 1)]
    X = [[_NEG] * (m + 1) for _ in range(n + 1)]   # gap in b (consumes a)
    Y = [[_NEG] * (m + 1) for _ in range(n + 1)]   # gap in a (consumes b)
    tM = [bytearray(m + 1) for _ in range(n + 1)]
    tX = [bytearray(m + 1) for _ in range(n + 1)]
    tY = [bytearray(m + 1) for _ in range(n + 1)]
    M[0][0] = 0.0
    for i in range(1, n + 1):
        X[i][0] = -(gap_open + (i - 1) * gap_extend)
        tX[i][0] = 0 if i == 1 else 1
    for j in range(1, m + 1):
        if semiglobal:
            M[0][j] = 0.0
        else:
            Y[0][j] = -(gap_open + (j - 1) * gap_extend)
            tY[0][j] = 0 if j == 1 else 2
    for i in range(1, n + 1):
        ai = a[i - 1]
        Mi, Xi, Yi = M[i], X[i], Y[i]
        Mp, Xp, Yp = M[i - 1], X[i - 1], Y[i - 1]
        tMi, tXi, tYi = tM[i], tX[i], tY[i]
        for j in range(1, m + 1):
            s = match if ai == b[j - 1] else mismatch
            pm, px, py = Mp[j - 1], Xp[j - 1], Yp[j - 1]
            if pm >= px and pm >= py:
                Mi[j], tMi[j] = pm + s, 0
            elif px >= py:
                Mi[j], tMi[j] = px + s, 1
            else:
                Mi[j], tMi[j] = py + s, 2
            o_m, e_x, o_y = Mp[j] - gap_open, Xp[j] - gap_extend, Yp[j] - gap_open
            if o_m >= e_x and o_m >= o_y:
                Xi[j], tXi[j] = o_m, 0
            elif e_x >= o_y:
                Xi[j], tXi[j] = e_x, 1
            else:
                Xi[j], tXi[j] = o_y, 2
            o_m, o_x, e_y = Mi[j - 1] - gap_open, Xi[j - 1] - gap_open, Yi[j - 1] - gap_extend
            if o_m >= e_y and o_m >= o_x:
                Yi[j], tYi[j] = o_m, 0
            elif e_y >= o_x:
                Yi[j], tYi[j] = e_y, 2
            else:
                Yi[j], tYi[j] = o_x, 1

    if semiglobal:
        best, end_j, state = _NEG, m, 0
        for j in range(1, m + 1):
            for st, mat in ((0, M), (1, X), (2, Y)):
                if mat[n][j] > best:
                    best, end_j, state = mat[n][j], j, st
        i, j = n, end_j
    else:
        i, j = n, m
        best, state = max((M[n][m], 0), (X[n][m], 1), (Y[n][m], 2))
    end_b = j
    out_a, out_b = [], []
    while i > 0 or (j > 0 and not semiglobal):
        if state == 0:
            prev = tM[i][j]
            out_a.append(a[i - 1]); out_b.append(b[j - 1])
            i -= 1; j -= 1
        elif state == 1:
            prev = tX[i][j]
            out_a.append(a[i - 1]); out_b.append("-")
            i -= 1
        else:
            prev = tY[i][j]
            out_a.append("-"); out_b.append(b[j - 1])
            j -= 1
        state = prev
    return Alignment("".join(reversed(out_a)), "".join(reversed(out_b)), best, j, end_b)


@dataclass
class PrimerSite:
    strand: str            # '+' primer matches top strand (forward), '-' bottom (reverse)
    start: int             # 1-based, top-strand coordinates
    end: int
    aln: Alignment
    mismatches_total: int
    mismatches_3p5: int    # mismatches/gaps among the 5 primer bases at the 3' end
    terminal_mismatch: bool


def _three_prime_errors(aln: Alignment, window: int = 5) -> tuple[int, bool]:
    seen = errors = 0
    terminal = False
    for x, y in zip(reversed(aln.a), reversed(aln.b)):
        bad = x != y
        if x != "-":
            if seen == 0 and bad:
                terminal = True
            seen += 1
        if bad:
            errors += 1
        if seen >= window:
            break
    return errors, terminal


def primer_binding(primer: str, template: str) -> list[PrimerSite]:
    """Best binding site of a primer on each template strand (sorted best first)."""
    n = len(template)
    res = []
    for strand, target in (("+", template), ("-", reverse_complement(template))):
        aln = align(primer, target, semiglobal=True)
        e3, term = _three_prime_errors(aln)
        if strand == "+":
            s, e = aln.b_start + 1, aln.b_end
        else:
            s, e = n - aln.b_end + 1, n - aln.b_start
        res.append(PrimerSite(strand, s, e, aln, aln.mismatches + aln.gaps, e3, term))
    res.sort(key=lambda r: r.aln.score, reverse=True)
    return res
