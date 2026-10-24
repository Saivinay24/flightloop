"""Read a 2D .su2 mesh and report quad signed areas (orientation) and counts."""
import numpy as np


def read(path):
    lines = open(path).read().split("\n")
    i = 0; elems = []; pts = []
    while i < len(lines):
        l = lines[i]
        if l.startswith("NELEM="):
            n = int(l.split("=")[1]); elems = [list(map(int, x.split()[1:5])) for x in lines[i + 1:i + 1 + n] if x.split()[0] == "9"]; i += n
        elif l.startswith("NPOIN="):
            n = int(l.split("=")[1].split()[0]); pts = np.array([list(map(float, x.split()[:2])) for x in lines[i + 1:i + 1 + n]]); i += n
        i += 1
    return np.array(elems), pts


def areas(path):
    E, P = read(path)
    a, b, c, d = (P[E[:, k]] for k in range(4))
    cr = lambda u, v: u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    return 0.5 * cr(c - a, d - b)
