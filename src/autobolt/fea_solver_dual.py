import tempfile

import build123d
import fenics
import gmsh
import meshio

import numpy as np

BOLT_TAG = 101
PLATE_TAG = 102

def import_step_and_capture_new_volumes(step_path: str) -> list[int]:
    """Import a STEP file and return the list of volume entity tags created by that import."""
    before = set(t for _, t in gmsh.model.getEntities(3))
    gmsh.merge(step_path)
    gmsh.model.occ.synchronize()
    after = set(t for _, t in gmsh.model.getEntities(3))
    new_tags = sorted(after - before)
    return new_tags

def _run_fenics_simulation(
        xdmf_file_mesh: str, 
        xdmf_file_surface_tags: str, 
        elastic_modulus: float,
        poissons_ratio: float,
        bolt_yield_strength: float,
        plate_yield_strength: float,
        traction_values: list[tuple],
        zero_disp_tags: list[int],
        traction_tags: list[int]
    ):
    mesh = fenics.cpp.mesh.Mesh()
    with fenics.cpp.io.XDMFFile(xdmf_file_mesh) as infile:
        infile.read(mesh)

    mvc_surf = fenics.MeshValueCollection("size_t", mesh, 2)
    with fenics.cpp.io.XDMFFile(xdmf_file_surface_tags) as infile:
        infile.read(mvc_surf, "name_to_read")
    mf = fenics.cpp.mesh.MeshFunctionSizet(mesh, mvc_surf)

    mvc_cell = fenics.MeshValueCollection("size_t", mesh, 3)
    with fenics.XDMFFile(xdmf_file_mesh) as infile:
        infile.read(mvc_cell, "name_to_read")
    cf = fenics.MeshFunction("size_t", mesh, mvc_cell)

    V = fenics.VectorFunctionSpace(mesh, "CG", 2)

    bcs = [
        fenics.DirichletBC(V, fenics.Constant((0.0, 0.0, 0.0)), mf, tag)
        for tag in zero_disp_tags
    ]

    mu = elastic_modulus / (2.0 * (1.0 + poissons_ratio))
    lmbda = elastic_modulus * poissons_ratio / (
        (1.0 + poissons_ratio) * (1.0 - 2.0 * poissons_ratio)
    )

    def epsilon(u):
        return 0.5 * (fenics.grad(u) + fenics.grad(u).T)

    def sigma(u):
        return 2.0 * mu * epsilon(u) + lmbda * fenics.tr(epsilon(u)) * fenics.Identity(3)
    
    u = fenics.TrialFunction(V)
    v = fenics.TestFunction(V)
    f = fenics.Constant((0.0, 0.0, 0.0))

    a = fenics.inner(sigma(u), epsilon(v)) * fenics.dx
    ds = fenics.Measure("ds", domain=mesh, subdomain_data=mf)
    L = fenics.dot(f, v) * fenics.dx

    for tag, traction in zip(traction_tags, traction_values):
        L += fenics.dot(fenics.Constant(traction), v) * ds(tag)

    u_sol = fenics.Function(V)
    fenics.solve(a == L, u_sol, bcs)

    s_dev = sigma(u_sol) - (1.0 / 3.0) * fenics.tr(sigma(u_sol)) * fenics.Identity(3)
    W0 = fenics.FunctionSpace(mesh, "DG", 0)
    vm0 = fenics.project(fenics.sqrt(3.0 / 2.0 * fenics.inner(s_dev, s_dev)), W0)
    vals = vm0.vector().get_local()
    markers = cf.array()

    bolt_vals = vals[markers == BOLT_TAG]
    plate_vals = vals[markers == PLATE_TAG]

    max_bolt_vm = float(np.max(bolt_vals))
    max_plate_vm = float(np.max(plate_vals))

    bolt_fos = float(bolt_yield_strength / max_bolt_vm)
    plate_fos = float(plate_yield_strength / max_plate_vm)

    return bolt_fos, plate_fos

def _calculate_fos_from_build123d_dual(
    plateA: build123d.Shape,
    plateB: build123d.Shape,
    bolt_cylinders: build123d.Shape,
    elastic_modulus: float,
    poissons_ratio: float,
    bolt_yield_strength: float,
    plate_yield_strength: float,
    traction_values: list[tuple],
):
    tempdir = tempfile.TemporaryDirectory()
    
    xdmf_file_mesh = tempdir.name + "/mesh.xdmf"
    xdmf_file_surface_tags = tempdir.name + "surface_tags.xdmf"
    mesh_file, zero_disp_tags, traction_tags = generate_tagged_mesh(plateA, plateB, bolt_cylinders, tempdir)

    convert_to_xdmf(mesh_file, "tetra", xdmf_file_mesh)
    convert_to_xdmf(mesh_file, "triangle", xdmf_file_surface_tags)

    return _run_fenics_simulation(
        xdmf_file_mesh,
        xdmf_file_surface_tags,
        elastic_modulus,
        poissons_ratio,
        bolt_yield_strength,
        plate_yield_strength,
        traction_values,
        zero_disp_tags,
        traction_tags
    )

def convert_to_xdmf(mesh_file, cell_type, out_path):
    mesh = meshio.read(mesh_file)

    cells = mesh.get_cells_type(cell_type)
    cell_data = mesh.get_cell_data("gmsh:physical", cell_type)
    out = meshio.Mesh(
        points=mesh.points,
        cells={cell_type: cells},
        cell_data={"name_to_read": [cell_data]},
    )
    meshio.write(out_path, out)

def generate_tagged_mesh(
    plateA: build123d.Shape,
    plateB: build123d.Shape,
    bolt_cylinders: build123d.Shape,
    tempdir: tempfile.TemporaryDirectory
):
    gmsh.initialize()

    step_plateA = tempdir.name + "/plateA.step"
    step_plateB = tempdir.name + "/plateB.step"
    step_bolts = tempdir.name + "/bolts.step"

    build123d.export_step(plateA, step_plateA)
    build123d.export_step(plateB, step_plateB)
    build123d.export_step(bolt_cylinders, step_bolts)

    plateA_vols = import_step_and_capture_new_volumes(step_plateA)
    plateB_vols = import_step_and_capture_new_volumes(step_plateB)
    bolt_vols = import_step_and_capture_new_volumes(step_bolts)

    in_dimtags = [(3, t) for t in (plateA_vols + plateB_vols + bolt_vols)]
    _, out_map = gmsh.model.occ.fragment(in_dimtags, [])
    gmsh.model.occ.synchronize()

    gmsh.model.occ.removeAllDuplicates()
    gmsh.model.occ.synchronize()

    bolt_set_pre = set(bolt_vols)

    bolt_frag_vols: set[int] = set()
    plate_frag_vols: set[int] = set()

    for (src_dim, src_tag), frags in zip(in_dimtags, out_map):
        frag_vols = [t for dim, t in frags if dim == 3]
        if src_tag in bolt_set_pre:
            bolt_frag_vols.update(frag_vols)
        else:
            plate_frag_vols.update(frag_vols)

    gmsh.model.addPhysicalGroup(3, sorted(bolt_frag_vols), tag=BOLT_TAG)
    gmsh.model.setPhysicalName(3, BOLT_TAG, "bolt")
    gmsh.model.addPhysicalGroup(3, sorted(plate_frag_vols), tag=PLATE_TAG)
    gmsh.model.setPhysicalName(3, PLATE_TAG, "plate")

    surfaces = gmsh.model.getEntities(2)

    for (sdim, stag) in surfaces:
        gmsh.model.addPhysicalGroup(sdim, [stag], tag=stag)

    surface_info = []
    for dim, s in gmsh.model.getEntities(2):
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.occ.getBoundingBox(dim, s)
        surface_info.append((s, ymin, ymax))

    global_ymin = min(ymin for _, ymin, _ in surface_info)
    global_ymax = max(ymax for _, _, ymax in surface_info)

    def y_face(target_y: float, tol: float = 1e-6) -> int:
        for tag, ymin, ymax in surface_info:
            if abs(ymax - ymin) < tol and abs(ymin - target_y) < tol:
                return tag
        raise RuntimeError(f"Could not find planar Y-face at y={target_y}")
    
    trac_pos_tag = y_face(global_ymin)
    trac_neg_tag = y_face(global_ymax)

    zero_disp_tags = [trac_pos_tag]
    traction_tags = [trac_neg_tag]

    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    gmsh.option.setNumber("Mesh.Optimize", 1)
    gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)

    gmsh.model.mesh.generate(3)

    gmsh.model.mesh.optimize("Netgen")

    gmsh.model.mesh.removeDuplicateNodes()

    msh_file = tempdir.name + "/mesh_merged.msh"

    gmsh.write(msh_file)

    gmsh.finalize()

    return msh_file, zero_disp_tags, traction_tags
