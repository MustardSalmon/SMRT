"""Input card for characteristic-line tracking on the Gmsh mesh."""

from __future__ import annotations

import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT))

import RayTracing as smrt
from Driver import PhysicsObj
from Driver.TrackFile import TrackFile
from RayTracing.GmshParser import GmshFileTranslator


scan_mode = "full"  # "full" or "one_sixth"; one-sixth requires a symmetric mesh
azim = 24
spacing = 0.03
num_proc = 4


if __name__ == "__main__":
    mesh_path = Path(__file__).with_name("HexPin.msh").resolve()
    if not mesh_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {mesh_path}. Run model.py first.")

    start = time.perf_counter()
    mesh = GmshFileTranslator(mesh_path, geom_type="hex").extractDataFromGmshFile()

    mesh_splitter = None
    if scan_mode == "one_sixth":
        mesh_splitter = smrt.MeshTreat.MeshSplitter(mesh)
        mesh_splitter.split("1/6")
    elif scan_mode != "full":
        raise ValueError("scan_mode must be 'full' or 'one_sixth'")

    track_factory = smrt.Track.TrackFactory(
        mesh,
        spacing,
        azim,
        num_proc=num_proc,
        mesh_spliter=mesh_splitter,
    )
    track_factory.tracking()

    rtm = PhysicsObj.RayTracingModule(mesh, track_factory, mesh_splitter)
    rtm_map = PhysicsObj.RTMMap([[rtm]], 1.0, 1)
    track_path = mesh_path.with_name(f"HexPin-{scan_mode}-{spacing}-{azim}.trk")
    TrackFile(str(track_path), "alpha").export([rtm], [rtm_map])

    print(f"Tracking time = {time.perf_counter() - start:.3f} sec")
    print(f"Track file written to {track_path}")
