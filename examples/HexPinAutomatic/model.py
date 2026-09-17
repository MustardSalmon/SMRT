# -*- coding: utf-8 -*-
# @File     : model.py
# @Author   : CaoWei
# @Colleges : Zhejiang University
# @Email    : 12245028@zju.edu.cn
# @Time     : 2026/9/17 16:22

from __future__ import annotations

import math
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

from GmshGeometryModeling import Material, PinCellObj, Plot, auxiliary


start = time.perf_counter()

uo2 = Material.Material("uo2")
uo2.setColor(1.0, 1.0, 0.82)

zr = Material.Material("zr")
zr.setColor(1.0, 0.72, 0.82)

co2 = Material.Material("co2")
co2.setColor(0.66, 0.85, 0.92)

mats = Material.Materials([uo2, zr, co2])
mats.create()

PinPitch = 1.08
radius1 = 0.4225
radius2 = 0.487571
edge_width = PinPitch / math.sqrt(3.0)

pin_cell = PinCellObj.HexRodLikePinCell(edge_width, [radius1, radius2], [uo2, zr, co2], mesh_sizes=[0.1, 0.2, 0.3],
                                        mesh_mode='automatic')
pin_cell.orientation = "x"
pin_cell.create()

mesh_path = Path(__file__).with_name("HexPin.msh")
auxiliary.save(mesh_path, "hex")
print(f"Gmsh mesh written to {mesh_path.resolve()}")

plot_path = Plot.plot_mesh(
    mesh_path,
    materials=mats,
    geom_type="hex",
    output_format="svg",        # png, svg, pdf
)
print(f"Mesh plot written to {plot_path}")

print(f"Geometry and meshing time = {time.perf_counter() - start:.3f} sec")
