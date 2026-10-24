"""Structured gmsh meshes (transfinite + recombine) for the textbook verification cases, written as .su2.

wedge  : 10 deg ramp starting at x = 0 after a flat run of 0.5, ramp length 1.5, top at y = 1.5.
diamond: symmetric 5% diamond (unit chord) in an 8-block H-mesh, box x in [-2, 3], y in [-3, 3].
Refinement ratio 2 per direction between levels; geometric progressions use g^(1/2) per halving, so every coarse
node is a node of the next level (nested, systematic).
"""
import math
import gmsh


def _line(p, q, n, prog=1.0):
    l = gmsh.model.geo.addLine(p, q)
    gmsh.model.geo.mesh.setTransfiniteCurve(l, n + 1, "Progression", prog)
    return l


def _surf(c1, c2, c3, c4):
    cl = gmsh.model.geo.addCurveLoop([c1, c2, c3, c4])
    s = gmsh.model.geo.addPlaneSurface([cl])
    gmsh.model.geo.mesh.setTransfiniteSurface(s)
    gmsh.model.geo.mesh.setRecombine(2, s)
    return s


def wedge(path, level, theta_deg=10.0):
    f = 2 ** level
    nx_up, nx_r, ny = 20 * f, 80 * f, 80 * f
    gy = 1.02 ** (1 / f)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("wedge")
    t = math.tan(math.radians(theta_deg)); L, H = 1.5, 1.5
    P = [gmsh.model.geo.addPoint(*xy, 0) for xy in [(-0.5, 0), (0, 0), (L, L * t), (L, H), (0, H), (-0.5, H)]]
    l_up = _line(P[0], P[1], nx_up); l_ramp = _line(P[1], P[2], nx_r)
    l_out = _line(P[2], P[3], ny, gy); l_top_r = _line(P[3], P[4], nx_r); l_top_u = _line(P[4], P[5], nx_up)
    l_in = _line(P[5], P[0], ny, 1 / gy); l_mid = _line(P[1], P[4], ny, gy)
    s1 = _surf(l_up, l_mid, l_top_u, l_in)
    s2 = _surf(l_ramp, l_out, l_top_r, -l_mid)
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(1, [l_up, l_ramp], name="wall")
    gmsh.model.addPhysicalGroup(1, [l_out, l_top_r, l_top_u, l_in], name="farfield")
    gmsh.model.addPhysicalGroup(2, [s1, s2], name="fluid")
    gmsh.model.mesh.generate(2)
    gmsh.write(str(path))
    n = gmsh.model.mesh.getElements(2)[1][0].size
    gmsh.finalize()
    return n


def diamond(path, level, tc=0.05):
    f = 2 ** level
    nu, nf, nw, ny = 24 * f, 40 * f, 40 * f, 60 * f
    gy = 1.06 ** (1 / f)
    gx_u = 1.05 ** (1 / f)
    gmsh.initialize(); gmsh.option.setNumber("General.Terminal", 0); gmsh.model.add("diamond")
    h = tc / 2
    pts = {"A": (-2, 0), "LE": (0, 0), "T": (0.5, h), "TE": (1, 0), "E": (3, 0), "Bm": (0.5, -h),
           "B": (-2, 3), "P": (0, 3), "C": (0.5, 3), "Q": (1, 3), "D": (3, 3),
           "Bl": (-2, -3), "Pl": (0, -3), "Cl": (0.5, -3), "Ql": (1, -3), "Dl": (3, -3)}
    G = {k: gmsh.model.geo.addPoint(x, y, 0) for k, (x, y) in pts.items()}
    L = {}
    def ln(a, b, n, prog=1.0):
        L[(a, b)] = _line(G[a], G[b], n, prog); return L[(a, b)]
    # horizontal-ish curves
    ln("A", "LE", nu, 1 / gx_u); ln("LE", "T", nf); ln("T", "TE", nf); ln("TE", "E", nw, gx_u)
    ln("LE", "Bm", nf); ln("Bm", "TE", nf)
    ln("B", "P", nu, 1 / gx_u); ln("P", "C", nf); ln("C", "Q", nf); ln("Q", "D", nw, gx_u)
    ln("Bl", "Pl", nu, 1 / gx_u); ln("Pl", "Cl", nf); ln("Cl", "Ql", nf); ln("Ql", "Dl", nw, gx_u)
    # vertical curves (cluster at the body / axis)
    for lo, hi in (("A", "B"), ("LE", "P"), ("T", "C"), ("TE", "Q"), ("E", "D")):
        ln(lo, hi, ny, gy)
    for lo, hi in (("A", "Bl"), ("LE", "Pl"), ("Bm", "Cl"), ("TE", "Ql"), ("E", "Dl")):
        ln(lo, hi, ny, gy)
    S = []
    def blk(b, r, t, l):  # bottom (left->right), right (up), top (left->right), left (up)
        S.append(_surf(L[b], L[r], -L[t], -L[l]))
    blk(("A", "LE"), ("LE", "P"), ("B", "P"), ("A", "B"))
    blk(("LE", "T"), ("T", "C"), ("P", "C"), ("LE", "P"))
    blk(("T", "TE"), ("TE", "Q"), ("C", "Q"), ("T", "C"))
    blk(("TE", "E"), ("E", "D"), ("Q", "D"), ("TE", "Q"))
    # lower half: 'up' direction is -y, keep orientation positive by reversing
    def blk_lo(top, r, bot, l):   # top = body-side curve (left->right), bot = outer (left->right), r/l from body outward
        cl = gmsh.model.geo.addCurveLoop([-L[top] if False else L[l], L[bot], -L[r], -L[top]])
        s = gmsh.model.geo.addPlaneSurface([cl]); gmsh.model.geo.mesh.setTransfiniteSurface(s); gmsh.model.geo.mesh.setRecombine(2, s); S.append(s)
    blk_lo(("A", "LE"), ("LE", "Pl"), ("Bl", "Pl"), ("A", "Bl"))
    blk_lo(("LE", "Bm"), ("Bm", "Cl"), ("Pl", "Cl"), ("LE", "Pl"))
    blk_lo(("Bm", "TE"), ("TE", "Ql"), ("Cl", "Ql"), ("Bm", "Cl"))
    blk_lo(("TE", "E"), ("E", "Dl"), ("Ql", "Dl"), ("TE", "Ql"))
    gmsh.model.geo.synchronize()
    gmsh.model.addPhysicalGroup(1, [L[("LE", "T")], L[("T", "TE")], L[("LE", "Bm")], L[("Bm", "TE")]], name="airfoil")
    far = [L[k] for k in [("B", "P"), ("P", "C"), ("C", "Q"), ("Q", "D"), ("Bl", "Pl"), ("Pl", "Cl"), ("Cl", "Ql"), ("Ql", "Dl"),
                          ("A", "B"), ("A", "Bl"), ("E", "D"), ("E", "Dl")]]
    gmsh.model.addPhysicalGroup(1, far, name="farfield")
    gmsh.model.addPhysicalGroup(2, S, name="fluid")
    gmsh.model.mesh.generate(2)
    gmsh.write(str(path))
    n = sum(e.size for e in gmsh.model.mesh.getElements(2)[1])
    gmsh.finalize()
    return n
