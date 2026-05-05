#!/usr/bin/env python3
"""Generate a minimal OpenFOAM polyMesh for the buflowRANS solver.

Creates a simple 2D-extruded structured mesh around a flat plate "airfoil"
with 5 boundaries: airfoil, empty, inlet, outlet, symmetry.
"""

import os, math

NX, NY = 6, 4
NZ = 1
DZ = 0.1

X_MIN, X_MAX = -1.0, 2.0
Y_MIN, Y_MAX = 0.0, 1.5

dx = (X_MAX - X_MIN) / NX
dy = (Y_MAX - Y_MIN) / NY

out_dir = "mesh/OFairfoilMesh"
os.makedirs(out_dir, exist_ok=True)

def pt_idx(ix, iy, iz):
    return iz * (NX + 1) * (NY + 1) + iy * (NX + 1) + ix

def cell_idx(ix, iy):
    return iy * NX + ix

nPts = (NX + 1) * (NY + 1) * (NZ + 1)
nCells = NX * NY * NZ

points = []
for iz in range(NZ + 1):
    for iy in range(NY + 1):
        for ix in range(NX + 1):
            x = X_MIN + ix * dx
            y = Y_MIN + iy * dy
            z = iz * DZ
            points.append((x, y, z))

internal_faces = []
owners_int = []
neighbours_int = []

for iy in range(NY):
    for ix in range(NX - 1):
        c0 = cell_idx(ix, iy)
        c1 = cell_idx(ix + 1, iy)
        p0 = pt_idx(ix + 1, iy, 0)
        p1 = pt_idx(ix + 1, iy + 1, 0)
        p2 = pt_idx(ix + 1, iy + 1, 1)
        p3 = pt_idx(ix + 1, iy, 1)
        internal_faces.append((p0, p3, p2, p1))
        owners_int.append(c0)
        neighbours_int.append(c1)

for iy in range(NY - 1):
    for ix in range(NX):
        c0 = cell_idx(ix, iy)
        c1 = cell_idx(ix, iy + 1)
        p0 = pt_idx(ix, iy + 1, 0)
        p1 = pt_idx(ix + 1, iy + 1, 0)
        p2 = pt_idx(ix + 1, iy + 1, 1)
        p3 = pt_idx(ix, iy + 1, 1)
        internal_faces.append((p0, p1, p2, p3))
        owners_int.append(c0)
        neighbours_int.append(c1)

nInternalFaces = len(internal_faces)

airfoil_faces, airfoil_owners = [], []
for ix in range(NX):
    c = cell_idx(ix, 0)
    p0 = pt_idx(ix, 0, 0)
    p1 = pt_idx(ix, 0, 1)
    p2 = pt_idx(ix + 1, 0, 1)
    p3 = pt_idx(ix + 1, 0, 0)
    airfoil_faces.append((p0, p1, p2, p3))
    airfoil_owners.append(c)

empty_faces_front, empty_owners_front = [], []
for iy in range(NY):
    for ix in range(NX):
        c = cell_idx(ix, iy)
        p0 = pt_idx(ix, iy, 1)
        p1 = pt_idx(ix + 1, iy, 1)
        p2 = pt_idx(ix + 1, iy + 1, 1)
        p3 = pt_idx(ix, iy + 1, 1)
        empty_faces_front.append((p0, p1, p2, p3))
        empty_owners_front.append(c)

empty_faces_back, empty_owners_back = [], []
for iy in range(NY):
    for ix in range(NX):
        c = cell_idx(ix, iy)
        p0 = pt_idx(ix, iy, 0)
        p1 = pt_idx(ix, iy + 1, 0)
        p2 = pt_idx(ix + 1, iy + 1, 0)
        p3 = pt_idx(ix + 1, 0 if iy == 0 else iy, 0)
        p0_ = pt_idx(ix, iy, 0)
        p1_ = pt_idx(ix, iy + 1, 0)
        p2_ = pt_idx(ix + 1, iy + 1, 0)
        p3_ = pt_idx(ix + 1, iy, 0)
        empty_faces_back.append((p0_, p1_, p2_, p3_))
        empty_owners_back.append(c)

empty_faces = empty_faces_front + empty_faces_back
empty_owners = empty_owners_front + empty_owners_back

inlet_faces, inlet_owners = [], []
for iy in range(NY):
    c = cell_idx(0, iy)
    p0 = pt_idx(0, iy, 0)
    p1 = pt_idx(0, iy + 1, 0)
    p2 = pt_idx(0, iy + 1, 1)
    p3 = pt_idx(0, iy, 1)
    inlet_faces.append((p0, p1, p2, p3))
    inlet_owners.append(c)

outlet_faces, outlet_owners = [], []
for iy in range(NY):
    c = cell_idx(NX - 1, iy)
    p0 = pt_idx(NX, iy, 0)
    p1 = pt_idx(NX, iy, 1)
    p2 = pt_idx(NX, iy + 1, 1)
    p3 = pt_idx(NX, iy + 1, 0)
    outlet_faces.append((p0, p1, p2, p3))
    outlet_owners.append(c)

symmetry_faces, symmetry_owners = [], []
for ix in range(NX):
    c = cell_idx(ix, NY - 1)
    p0 = pt_idx(ix, NY, 0)
    p1 = pt_idx(ix + 1, NY, 0)
    p2 = pt_idx(ix + 1, NY, 1)
    p3 = pt_idx(ix, NY, 1)
    symmetry_faces.append((p0, p1, p2, p3))
    symmetry_owners.append(c)

all_faces = internal_faces
all_owners = owners_int[:]

bdry_start = {}
bdry_nfaces = {}

bdry_start["airfoil"] = len(all_faces)
all_faces += airfoil_faces
all_owners += airfoil_owners
bdry_nfaces["airfoil"] = len(airfoil_faces)

bdry_start["empty"] = len(all_faces)
all_faces += empty_faces
all_owners += empty_owners
bdry_nfaces["empty"] = len(empty_faces)

bdry_start["inlet"] = len(all_faces)
all_faces += inlet_faces
all_owners += inlet_owners
bdry_nfaces["inlet"] = len(inlet_faces)

bdry_start["outlet"] = len(all_faces)
all_faces += outlet_faces
all_owners += outlet_owners
bdry_nfaces["outlet"] = len(outlet_faces)

bdry_start["symmetry"] = len(all_faces)
all_faces += symmetry_faces
all_owners += symmetry_owners
bdry_nfaces["symmetry"] = len(symmetry_faces)

nFaces = len(all_faces)

with open(os.path.join(out_dir, "points"), "w") as f:
    f.write(f"{nPts}\n(\n")
    for p in points:
        f.write(f"({p[0]} {p[1]} {p[2]})\n")
    f.write(")\n")

with open(os.path.join(out_dir, "faces"), "w") as f:
    f.write(f"{nFaces}\n(\n")
    for face in all_faces:
        f.write(f"{len(face)}({' '.join(str(v) for v in face)})\n")
    f.write(")\n")

with open(os.path.join(out_dir, "owner"), "w") as f:
    f.write(f"{nFaces}\n(\n")
    for o in all_owners:
        f.write(f"{o}\n")
    f.write(")\n")

with open(os.path.join(out_dir, "neighbour"), "w") as f:
    f.write(f"{nInternalFaces}\n(\n")
    for n in neighbours_int:
        f.write(f"{n}\n")
    f.write(")\n")

boundary_order = ["airfoil", "empty", "inlet", "outlet", "symmetry"]
with open(os.path.join(out_dir, "boundary"), "w") as f:
    f.write(f"{len(boundary_order)}\n(\n")
    for name in boundary_order:
        f.write(f"    {name}\n")
        f.write("    {\n")
        f.write(f"        type patch\n")
        f.write(f"        nFaces {bdry_nfaces[name]}\n")
        f.write(f"        startFace {bdry_start[name]}\n")
        f.write("    }\n")
    f.write(")\n")

print(f"Mesh generated: {nCells} cells, {nFaces} faces, {nPts} points")
print(f"Internal faces: {nInternalFaces}")
for name in boundary_order:
    print(f"  {name}: nFaces={bdry_nfaces[name]}, startFace={bdry_start[name]}")
