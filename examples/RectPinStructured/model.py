# -*- coding: utf-8 -*-
# @File     : model.py
# @Author   : CaoWei
# @Colleges : Zhejiang University
# @Email    : 12245028@zju.edu.cn
# @Time     : 2026/9/17 16:22

from __future__ import annotations

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

h2o = Material.Material("h2o")
h2o.setColor(0.66, 0.85, 0.92)

mats = Material.Materials([uo2, zr, h2o])
mats.create()

PinPitch = 1.26
radius1 = 0.46
radius2 = 0.54

pin_cell = PinCellObj.RectRodLikePinCell(PinPitch, PinPitch, [radius1, radius2], [uo2, zr, h2o],
                                         radial_meshes=[3, 1, 2], azimuthal_meshes=[24, 24, 24], mesh_mode='structured')
pin_cell.create()

mesh_path = Path(__file__).with_name("RectPin.msh")
auxiliary.save(mesh_path, "rect")
print(f"Gmsh mesh written to {mesh_path.resolve()}")

plot_path = Plot.plot_mesh(
    mesh_path,
    materials=mats,
    geom_type="rect",
    output_format="svg",        # png, svg, pdf
)
print(f"Mesh plot written to {plot_path}")

print(f"Geometry and meshing time = {time.perf_counter() - start:.3f} sec")
