import build123d


def create_two_plate_assembly_components(
    plate_length_m: float,
    plate_width_m: float,
    plate_thickness_m: float,
    num_holes: int,
    hole_radius_m: float,
    edge_margin_m: float,
    hole_spacing_m: float,
    hole_offset_from_bottom_m: float,
    plate_gap_mm: float,
):
    """
    Create a two-plate assembly with through-holes and bolt-like cylinders.

    Parameters:
        plate_length_m (float): Length of the plates in meters.
        plate_width_m (float): Width of the plates in meters.
        plate_thickness_m (float): Thickness of the plates in meters.
        num_holes (int): Number of through-holes in each plate.
        hole_radius_m (float): Radius of each through-hole in meters.
        edge_margin_m (float): Margin from the edge of the plate to the first hole in meters.
        hole_spacing_m (float): Spacing between the centers of adjacent holes in meters.
        hole_offset_from_bottom_m (float): Vertical offset of the holes from the bottom edge of the plate in meters.
        plate_gap_mm (float): Gap between the two plates in millimeters

    Returns:
        build123d.Compound: The assembled two-plate structure with bolt-like cylinders.

    """
    # convert the gap to meters
    gap_m = plate_gap_mm * 1e-4

    # 1) Build Plate A (with through‑holes)
    with build123d.BuildPart() as bp:
        build123d.Box(plate_length_m, plate_width_m, plate_thickness_m)
        # center the box around the XY origin
        bp.part = bp.part.move(build123d.Pos(plate_length_m / 2, plate_width_m / 2, 0))

        # compute X‑coordinates of hole centers
        x0 = plate_length_m - edge_margin_m
        x_coords = [x0 - i * hole_spacing_m for i in range(num_holes)]

        # subtract each hole
        for x in x_coords:
            with build123d.Locations(build123d.Pos(x, hole_offset_from_bottom_m, 0)):
                build123d.Cylinder(
                    radius=hole_radius_m,
                    height=plate_thickness_m + 1e-3,  # +1 mm to ensure clean cut
                    mode=build123d.Mode.SUBTRACT,
                )
    plateA = bp.part

    # 2) Copy Plate A → shift up by thickness+gap → mirror in Y → re‑position
    plateB = build123d.copy(plateA)
    # shift in Z by plate_thickness + gap
    plateB = plateB.locate(
        build123d.Location((0, 0, plate_thickness_m + gap_m), (0, 0, 1), 0)
    )
    # mirror about the XZ plane (so y→−y)
    plateB = build123d.mirror(plateB, about=build123d.Plane.XZ)
    # then translate back up in Y by 2*hole_center_y to line holes up
    plateB = plateB.locate(
        build123d.Location((0, 2 * hole_offset_from_bottom_m, 0), (0, 0, 1), 0)
    )

    # 3) Create the bolt‑like cylinders (conformal, centered in the holes)
    bolt_height_m = 2 * plate_thickness_m + gap_m
    z_offset_m = plate_thickness_m / 2 + gap_m / 2  # start mid‑thickness of Plate A

    with build123d.BuildPart() as bp2:
        for x in x_coords:
            with build123d.Locations(
                build123d.Pos(x, hole_offset_from_bottom_m, z_offset_m)
            ):
                build123d.Cylinder(
                    radius=hole_radius_m,
                    height=bolt_height_m,
                    align=(
                        build123d.Align.CENTER,
                        build123d.Align.CENTER,
                        build123d.Align.CENTER,
                    ),
                )
    bolt_cylinders = bp2.part

    # 4) Assemble and return
    return plateA, plateB, bolt_cylinders


def create_two_plate_assembly(
    plate_length_m: float,
    plate_width_m: float,
    plate_thickness_m: float,
    num_holes: int,
    hole_radius_m: float,
    edge_margin_m: float,
    hole_spacing_m: float,
    hole_offset_from_bottom_m: float,
    plate_gap_mm: float,
):
    plateA, plateB, bolt_cylinders = create_two_plate_assembly_components(
        plate_length_m,
        plate_width_m,
        plate_thickness_m,
        num_holes,
        hole_radius_m,
        edge_margin_m,
        hole_spacing_m,
        hole_offset_from_bottom_m,
        plate_gap_mm,
    )

    assembly = build123d.Compound(children=[plateA, plateB, bolt_cylinders])
    return assembly