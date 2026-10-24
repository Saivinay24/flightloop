"""Minimal reader for SU2's raw-appended binary VTU files (header_type UInt64)."""
import re
import numpy as np

DT = {"Float32": np.float32, "Float64": np.float64, "Int32": np.int32, "UInt8": np.uint8, "Int64": np.int64}


def read_vtu(path):
    raw = open(path, "rb").read()
    i = raw.index(b"<AppendedData")
    j = raw.index(b"_", i) + 1
    head = raw[:i].decode("latin1")
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents= ?"(\d+)" offset="(\d+)"', head):
        typ, name, nc, off = m.group(1), m.group(2) or "Points", int(m.group(3)), int(m.group(4))
        n = int(np.frombuffer(raw[j + off: j + off + 8], np.uint64)[0])
        arr = np.frombuffer(raw[j + off + 8: j + off + 8 + n], DT[typ])
        out[name] = arr.reshape(-1, nc) if nc > 1 else arr
    return out
