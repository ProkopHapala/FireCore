#!/usr/bin/env python3
"""Fail-loud checks for XH-neighbor alignment (terrace vs edge). Wulff COM support is not replaced."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT))

from pyBall.FFfit_utils import (
    facet_kind_from_vec, halogen_symbol_for_nbhd_tag, heavy_neighbor_lists, heavies_near_110_extrema,
    is_xh_align_terrace, is_xh_on_miller_111, miller_110_unit_normals, miller_110_unsigned_axes,
    miller_111_unit_normals, neighborhood_xh_groups, write_hydride_halogen_xyz,
    xh_bonds_from_coords, xh_bonds_from_topology, xh_miller_111_cosine, xh_unit_dirs,
)

DFTB_OCTA = Path("/home/prokop/SIMULATIONS/SiNCs/DFTB/octahedron_C/geo_end.xyz")
A_CC = 1.54
R_CH = 1.09


def _xyz_Zpos(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    smap = {"H": 1, "C": 6, "Si": 14}
    pos = np.zeros((n, 3)); Z = np.zeros(n, dtype=np.int32)
    for i in range(n):
        w = lines[2 + i].split()
        Z[i] = smap[w[0]]
        pos[i] = (float(w[1]), float(w[2]), float(w[3]))
    return pos, Z


def _chain(h2):
    """C0–C1–C2, one H each. h2 = H on C2 relative to (2a, 0, R_CH)."""
    pos = np.array([
        [0.0, 0.0, 0.0], [A_CC, 0.0, 0.0], [2.0 * A_CC, 0.0, 0.0],
        [0.0, 0.0, R_CH], [A_CC, 0.0, R_CH], np.array([2.0 * A_CC, 0.0, R_CH]) + np.asarray(h2, dtype=np.float64),
    ])
    Z = np.array([6, 6, 6, 1, 1, 1], dtype=np.int32)
    bonds = np.array([[0, 1], [1, 2], [0, 3], [1, 4], [2, 5]], dtype=np.int32)
    xh, nH = xh_bonds_from_topology(Z, bonds)
    adj = heavy_neighbor_lists(pos, Z, bonds_ij=bonds)
    dirs = xh_unit_dirs(pos, xh)
    return pos, Z, xh, nH, adj, dirs, bonds


def test_facet_kind_from_vec_still_wulff():
    """Revert path: {111} octa families still score a [111] vector as 111, a [110] as 110."""
    assert facet_kind_from_vec([1.0, 1.0, 1.0], families=('111', '110')) == '111'
    assert facet_kind_from_vec([1.0, 1.0, 0.0], families=('111', '110')) == '110'
    assert facet_kind_from_vec([1.0, 0.0, 0.0], families=('100',)) == '100'


def test_parallel_ch_chain_all_terrace():
    pos, Z, xh, nH, adj, dirs, bonds = _chain([0.0, 0.0, 0.0])
    for i in (0, 1, 2):
        if not is_xh_align_terrace(i, nH, adj, dirs, align_cos=0.9, max_misaligned=1):
            raise AssertionError(f"parallel CH: heavy {i} should be terrace")
    nbhd, _ = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=bonds, families=('111', '110'), facet_mode='xh_align', face_families=('111',))
    if any(k.endswith('@edge') and len(nbhd[k]) for k in nbhd):
        raise AssertionError(f"parallel chain should have no @edge, got { {k: len(v) for k,v in nbhd.items() if len(v)} }")


def test_kinked_end_is_edge_middle_stays_terrace():
    """End CH bent ~50° from the group: that site is the ridge. Middle has one bad neighbor → still terrace."""
    pos, Z, xh, nH, adj, dirs, bonds = _chain([0.80, 0.0, -0.40])
    if is_xh_align_terrace(2, nH, adj, dirs, align_cos=0.9, max_misaligned=1):
        raise AssertionError("kinked end C2 should be edge")
    if not is_xh_align_terrace(1, nH, adj, dirs, align_cos=0.9, max_misaligned=1):
        raise AssertionError("middle C1 with one misaligned neighbor should stay terrace")
    if not is_xh_align_terrace(0, nH, adj, dirs, align_cos=0.9, max_misaligned=1):
        raise AssertionError("far end C0 (aligned to middle) should stay terrace")
    nbhd, _ = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=bonds, families=('111', '110'), facet_mode='xh_align', face_families=('111',))
    edge = {k: np.asarray(v) for k, v in nbhd.items() if k.endswith('@edge') and len(v)}
    heavies_edge = set()
    for v in edge.values():
        heavies_edge.update(int(i) for i, _j in v)
    if 2 not in heavies_edge:
        raise AssertionError(f"C2 missing from @edge tags { {k: len(v) for k,v in nbhd.items() if len(v)} }")
    if 1 in heavies_edge:
        raise AssertionError("C1 (one row in) must not be tagged @edge")


def test_isolated_ch_is_terrace():
    """Diamond {111}: terrace CH has no hydride heavy-neighbors (only bulk C)."""
    pos = np.array([[0.0, 0.0, 0.0], [A_CC, 0.0, 0.0], [0.0, 0.0, R_CH]])
    Z = np.array([6, 6, 1], dtype=np.int32)
    bonds = np.array([[0, 1], [0, 2]], dtype=np.int32)
    xh, nH = xh_bonds_from_topology(Z, bonds)
    adj = heavy_neighbor_lists(pos, Z, bonds_ij=bonds)
    dirs = xh_unit_dirs(pos, xh)
    if int(nH[0]) != 1 or int(nH[1]) != 0:
        raise AssertionError(f"want nH=(1,0) got {nH[:2]}")
    if not is_xh_align_terrace(0, nH, adj, dirs):
        raise AssertionError("isolated CH (bulk neighbors only) must be terrace")


def test_halogen_symbol_face_vs_edge():
    """Jmol map: F/Cl = CH face/edge, Br/I = CH₂ face/edge. Miller face only if in face_families."""
    ff111 = ('111',)
    if halogen_symbol_for_nbhd_tag('XH@111', ff111) != 'F':
        raise AssertionError("CH@{111} face must be F")
    if halogen_symbol_for_nbhd_tag('XH@edge', ff111) != 'Cl':
        raise AssertionError("CH edge must be Cl")
    if halogen_symbol_for_nbhd_tag('XH2@edge', ff111) != 'I':
        raise AssertionError("CH₂ edge must be I")
    if halogen_symbol_for_nbhd_tag('XH@110', ff111) != 'Cl':
        raise AssertionError("CH@{110} is not an octa face → Cl (edge), not F")
    ff110 = ('110',)
    if halogen_symbol_for_nbhd_tag('XH@110', ff110) != 'F':
        raise AssertionError("CH@{110} on rhombic is the face → F")
    if halogen_symbol_for_nbhd_tag('XH2@110', ff110) != 'Br':
        raise AssertionError("CH₂@{110} face would be Br")
    if halogen_symbol_for_nbhd_tag('XH2@edge', ff110) != 'I':
        raise AssertionError("CH₂ edge on rhombic must be I")


def test_write_hydride_halogen_xyz(tmp_path):
    pos = np.zeros((2, 3)); pos[1, 2] = R_CH
    Z = np.array([6, 1], dtype=np.int32)
    p = tmp_path / "h.xyz"
    write_hydride_halogen_xyz(p, pos, Z, {0: 'XH@111', 1: 'XH@111'}, ('111',), comment='octa')
    lines = p.read_text().splitlines()
    if lines[0].strip() != '2' or not lines[2].startswith('C ') or not lines[3].startswith('F '):
        raise AssertionError(p.read_text())
    if 'F=XH-face' not in lines[1] or 'F:XH@111=1' not in lines[1]:
        raise AssertionError(lines[1])


def test_mixed_xh_vs_xh2_rim():
    """XH next to XH₂ is a rim. XH₂ with a single CH contact is a small {100} terrace (diamond has no CH₂–CH₂ bonds)."""
    pos = np.array([
        [0.0, 0.0, 0.0], [A_CC, 0.0, 0.0],
        [0.0, 0.0, R_CH], [A_CC, 0.0, R_CH], [A_CC, 0.0, -R_CH],
    ])
    Z = np.array([6, 6, 1, 1, 1], dtype=np.int32)
    bonds = np.array([[0, 1], [0, 2], [1, 3], [1, 4]], dtype=np.int32)
    xh, nH = xh_bonds_from_topology(Z, bonds)
    adj = heavy_neighbor_lists(pos, Z, bonds_ij=bonds)
    dirs = xh_unit_dirs(pos, xh)
    if int(nH[0]) != 1 or int(nH[1]) != 2:
        raise AssertionError(f"want nH=(1,2) got {nH[:2]}")
    if is_xh_align_terrace(0, nH, adj, dirs):
        raise AssertionError("XH with only an XH₂ neighbor is a rim, not a {111} terrace")
    if not is_xh_align_terrace(1, nH, adj, dirs, xh2_rim_terrace=True):
        raise AssertionError("XH₂ with one CH rim contact must stay a {100} terrace")
    if is_xh_align_terrace(1, nH, adj, dirs, xh2_rim_terrace=False):
        raise AssertionError("without {100} rim rule, CH₂ next to CH is a vertex, not a face")


def _nbhd_n(nbhd):
    return {k: int(len(v)) for k, v in nbhd.items() if len(v)}


@pytest.mark.skipif(not DFTB_OCTA.is_file(), reason="L2 octa DFTB geo_end.xyz not on this disk")
def test_l2_octa_xh_align_vs_wulff_counts():
    """Same hydrides; xh_align must not use Wulff {110} ribbon as the face/edge split. Print both counts."""
    pos, Z = _xyz_Zpos(DFTB_OCTA)
    xh, nH = xh_bonds_from_coords(pos, Z)
    n_xh = int(xh['XH'].shape[0])
    n_xh2 = int(xh['XH2'].shape[0])
    kw = dict(pos=pos, Z=Z, xh_groups=xh, nH=nH, bonds_ij=None, families=('111', '110'))
    wulff, _ = neighborhood_xh_groups(**kw, facet_mode='wulff')
    align, _ = neighborhood_xh_groups(**kw, facet_mode='xh_align', face_families=('111',))
    def _xh_sum(nb, cls):
        return sum(int(len(v)) for k, v in nb.items() if k.startswith(cls + '@'))
    if _xh_sum(wulff, 'XH') != n_xh or _xh_sum(align, 'XH') != n_xh:
        raise AssertionError(f"XH bond partition lost: n={n_xh} wulff={_xh_sum(wulff,'XH')} align={_xh_sum(align,'XH')}")
    if _xh_sum(wulff, 'XH2') != n_xh2 or _xh_sum(align, 'XH2') != n_xh2:
        raise AssertionError(f"XH2 bond partition lost: n={n_xh2} wulff={_xh_sum(wulff,'XH2')} align={_xh_sum(align,'XH2')}")
    n111_w = int(len(wulff.get('XH@111', ())))
    n110_w = int(len(wulff.get('XH@110', ())))
    n111_a = int(len(align.get('XH@111', ())))
    nedge_a = int(len(align.get('XH@edge', ())))
    print(f"L2 octahedron_C  nH={int((Z==1).sum())}  XH={n_xh} XH2={n_xh2}")
    print(f"  wulff     {_nbhd_n(wulff)}")
    print(f"  xh_align  {_nbhd_n(align)}")
    if n110_w < 1:
        raise AssertionError("Wulff octa should still tag some CH@{110}; do not delete facet_kind_from_vec")
    if 'XH@110' in align and len(align['XH@110']):
        raise AssertionError("xh_align on octa (face_families={111}) must tag non-terrace as @edge, not leftover @110")
    if n111_a < 1:
        raise AssertionError(f"xh_align produced no CH@{{111}} terrace (n111={n111_a} nedge={nedge_a})")
    if n111_a <= n111_w:
        raise AssertionError(f"xh_align should keep the Wulff {{110}} ribbon as {{111}} (align {n111_a} vs wulff {n111_w})")
    if nedge_a < 1:
        raise AssertionError("xh_align should still mark the geometric ridge as CH edge")
    adj = heavy_neighbor_lists(pos, Z)
    n_iso = 0
    for i in range(len(Z)):
        if int(nH[i]) != 1:
            continue
        if not any(int(nH[j]) > 0 for j in adj[i]):
            n_iso += 1
    if n111_a != n_iso:
        raise AssertionError(f"octa CH@{{111}} should be the isolated CH (no hydride neighbors): {n111_a} vs iso {n_iso}")


DFTB_RHOMB = Path("/home/prokop/SIMULATIONS/SiNCs/DFTB/rhombic_dodecahedron_C/geo_end.xyz")

@pytest.mark.skipif(not DFTB_RHOMB.is_file(), reason="L2 rhombic DFTB geo_end.xyz not on this disk")
def test_l2_rhombic_ch2_is_edge_not_110_face():
    """{110} is the rhombic face (CH). CH₂ are vertices — not CH₂@{110} face."""
    pos, Z = _xyz_Zpos(DFTB_RHOMB)
    xh, nH = xh_bonds_from_coords(pos, Z)
    align, _ = neighborhood_xh_groups(pos, Z, xh, nH, bonds_ij=None, families=('110',), facet_mode='xh_align', face_families=('110',))
    print(f"L2 rhombic  xh_align {_nbhd_n(align)}")
    if align.get('XH2@110') is not None and len(align['XH2@110']):
        raise AssertionError(f"rhombic CH₂ tagged as {{110}} face: {len(align['XH2@110'])} bonds")
    n_xh2_edge = int(len(align.get('XH2@edge', ())))
    if n_xh2_edge != int(xh['XH2'].shape[0]):
        raise AssertionError(f"all rhombic CH₂ bonds should be @edge, got {n_xh2_edge} vs {xh['XH2'].shape[0]}")
    if int(len(align.get('XH@110', ()))) < 1:
        raise AssertionError("rhombic should still have some isolated CH as {{110}} face")


def test_miller_111_cosine_gate():
    """Sitting face from COM octant; X–H must match that ⟨111⟩. ⟨110⟩/⟨100⟩ X–H on a {111} site are edge."""
    N = miller_111_unit_normals()
    if N.shape != (8, 3):
        raise AssertionError(N.shape)
    if abs(xh_miller_111_cosine([1.0, 1.0, 1.0]) - 1.0) > 1e-12:
        raise AssertionError(xh_miller_111_cosine([1.0, 1.0, 1.0]))
    c110 = xh_miller_111_cosine([1.0, 1.0, 0.0])
    if abs(c110 - np.sqrt(2.0 / 3.0)) > 1e-9:
        raise AssertionError(c110)
    r111 = [1.0, 1.0, 1.0]
    if not is_xh_on_miller_111([1.0, 1.0, 1.0], r111, align_cos=0.90):
        raise AssertionError("X–H along the sitting ⟨111⟩ must be a face")
    if is_xh_on_miller_111([1.0, 1.0, 0.0], r111, align_cos=0.90):
        raise AssertionError("X–H along [110] on a {111} site must be edge at 0.90")
    if is_xh_on_miller_111([1.0, 0.0, 0.0], r111, align_cos=0.90):
        raise AssertionError("X–H along [100] on a {111} site must be edge")
    if is_xh_on_miller_111([-1.0, -1.0, -1.0], r111, align_cos=0.90):
        raise AssertionError("inward X–H on a {111} site must be edge")


@pytest.mark.skipif(not DFTB_OCTA.is_file() or not DFTB_RHOMB.is_file(), reason="L2 DFTB geo_end.xyz not on this disk")
def test_l2_miller_111_enlarges_faces_vs_xh_align():
    """X–H vs 8 ⟨111⟩ should tag more face CH than isolated-CH xh_align, especially rhombic."""
    def _xh_sum(nb, cls):
        return sum(int(len(v)) for k, v in nb.items() if k.startswith(cls + '@'))
    for label, path in (("octahedron_C", DFTB_OCTA), ("rhombic_dodecahedron_C", DFTB_RHOMB)):
        pos, Z = _xyz_Zpos(path)
        xh, nH = xh_bonds_from_coords(pos, Z)
        kw = dict(pos=pos, Z=Z, xh_groups=xh, nH=nH, bonds_ij=None, families=('111', '110'))
        align, _ = neighborhood_xh_groups(**kw, facet_mode='xh_align', face_families=('111',) if 'octa' in label else ('110',))
        mill, _ = neighborhood_xh_groups(**kw, facet_mode='miller_111', align_cos=0.90)
        n_xh = int(xh['XH'].shape[0])
        if _xh_sum(mill, 'XH') != n_xh:
            raise AssertionError(f"{label}: miller_111 lost XH bonds {_xh_sum(mill,'XH')} vs {n_xh}")
        n111_m = int(len(mill.get('XH@111', ())))
        n_iso_face = int(len(align.get('XH@111', ()))) if 'octa' in label else int(len(align.get('XH@110', ())))
        print(f"L2 {label}  miller_111 {_nbhd_n(mill)}")
        print(f"           xh_align   {_nbhd_n(align)}")
        if n111_m < 1:
            raise AssertionError(f"{label}: miller_111 produced no CH@{{111}}")
        if n111_m <= n_iso_face:
            raise AssertionError(f"{label}: miller_111 CH@{{111}}={n111_m} should exceed xh_align face {n_iso_face}")


def test_miller_110_axes_12_or_6():
    n12 = miller_110_unit_normals()
    n6 = miller_110_unsigned_axes()
    if n12.shape != (12, 3) or n6.shape != (6, 3):
        raise AssertionError(f"want (12,3) and (6,3) got {n12.shape} {n6.shape}")
    if np.max(np.abs(np.linalg.norm(n12, axis=1) - 1.0)) > 1e-12:
        raise AssertionError("signed ⟨110⟩ must be unit")
    if np.max(np.abs(np.linalg.norm(n6, axis=1) - 1.0)) > 1e-12:
        raise AssertionError("unsigned ⟨110⟩ must be unit")


@pytest.mark.skipif(not DFTB_OCTA.is_file(), reason="L2 octa DFTB geo_end.xyz not on this disk")
def test_heavies_near_110_extrema_slab():
    """Interior C (nearest COM) is not a ⟨110⟩ cap; some surface C are. 0.5 Å slab."""
    pos, Z = _xyz_Zpos(DFTB_OCTA)
    heav = np.flatnonzero(Z > 1)
    com = pos[heav].mean(axis=0)
    i = int(heav[np.argmin(np.linalg.norm(pos[heav] - com, axis=1))])
    near = heavies_near_110_extrema(pos, Z, below_A=0.5)
    if near[i]:
        raise AssertionError(f"COM heavy {i} tagged as ⟨110⟩ extremity")
    if not np.any(near[heav]):
        raise AssertionError("no C/Si at ⟨110⟩ extrema")
    n_cap = int(np.sum(near[heav]))
    print(f"L2 octa  ⟨110⟩-cap heavies {n_cap}/{heav.size}  COM atom {i} near={bool(near[i])}")


@pytest.mark.skipif(not DFTB_OCTA.is_file() or not DFTB_RHOMB.is_file(), reason="L2 DFTB geo_end.xyz not on this disk")
def test_l2_ridge_110_partitions_xh():
    """ridge_110 must keep every X–H bond. Print counts vs miller_111 (trial, no face-count claim)."""
    def _xh_sum(nb, cls):
        return sum(int(len(v)) for k, v in nb.items() if k.startswith(cls + '@'))
    for label, path, ff in (("octahedron_C", DFTB_OCTA, ('111',)), ("rhombic_dodecahedron_C", DFTB_RHOMB, ('110',))):
        pos, Z = _xyz_Zpos(path)
        xh, nH = xh_bonds_from_coords(pos, Z)
        kw = dict(pos=pos, Z=Z, xh_groups=xh, nH=nH, bonds_ij=None, families=('111', '110'), face_families=ff)
        ridge, _ = neighborhood_xh_groups(**kw, facet_mode='ridge_110', ridge_below_A=0.5)
        mill, _ = neighborhood_xh_groups(**kw, facet_mode='miller_111', align_cos=0.90)
        n_xh = int(xh['XH'].shape[0])
        if _xh_sum(ridge, 'XH') != n_xh:
            raise AssertionError(f"{label}: ridge_110 lost XH {_xh_sum(ridge,'XH')} vs {n_xh}")
        print(f"L2 {label}  ridge_110 {_nbhd_n(ridge)}")
        print(f"           miller_111 {_nbhd_n(mill)}")


