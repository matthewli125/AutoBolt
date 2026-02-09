from .fea_solver import _calculate_fos_from_build123d
from .parametric_cad_solver import create_two_plate_assembly, create_two_plate_assembly_components

from .fea_solver_dual import _calculate_fos_from_build123d_dual


def calculate_fos_dual(
    plate_length_m: float,
    plate_width_m: float,
    plate_thickness_m: float,
    num_holes: int,
    hole_radius_m: float,
    edge_margin_m: float,
    hole_spacing_m: float,
    hole_offset_from_bottom_m: float,
    plate_gap_mm: float,
    elastic_modulus: float,
    poissons_ratio: float,
    bolt_yield_strength: float,
    plate_yield_strength: float,
    traction_values: list[tuple],
) -> tuple[float, float]:
    
    plateA, plateB, bolt_cylinders = create_two_plate_assembly(
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

    bolt_fos, plate_fos = _calculate_fos_from_build123d_dual(
        plateA,
        plateB,
        bolt_cylinders,
        elastic_modulus,
        poissons_ratio,
        bolt_yield_strength,
        plate_yield_strength,
        traction_values,
    )

    return bolt_fos, plate_fos

def calculate_fos(
    plate_length_m: float,
    plate_width_m: float,
    plate_thickness_m: float,
    num_holes: int,
    hole_radius_m: float,
    edge_margin_m: float,
    hole_spacing_m: float,
    hole_offset_from_bottom_m: float,
    plate_gap_mm: float,
    elastic_modulus: float,
    poissons_ratio: float,
    yield_strength: float,
    traction_values: list[tuple],
) -> float:
    """
    Calculate the Factor of Safety (FOS) for a two-plate assembly with bolt holes.

    Parameters:
        plate_length_m (float): Length of the plates in meters.
        plate_width_m (float): Width of the plates in meters.
        plate_thickness_m (float): Thickness of the plates in meters.
        num_holes (int): Number of bolt holes.
        hole_radius_m (float): Radius of each bolt hole in meters.
        edge_margin_m (float): Margin from the edge to the first hole in meters.
        hole_spacing_m (float): Spacing between the centers of adjacent holes in meters.
        hole_offset_from_bottom_m (float): Offset of the holes from the bottom edge in meters
        plate_gap_mm (float): Gap between the two plates in millimeters.
        elastic_modulus (float): Young's modulus of the material.
        poissons_ratio (float): Poisson's ratio of the material.
        yield_strength (float): Yield strength of the material.
        traction_values (list[tuple]): List of traction force vectors to apply on specified surfaces.

    Returns:
        float: The calculated Factor of Safety (FOS).
    """
    # Create the two-plate assembly
    assembly = create_two_plate_assembly(
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

    # Calculate the factor of safety
    fos = _calculate_fos_from_build123d(
        build123d_object=assembly,
        elastic_modulus=elastic_modulus,
        poissons_ratio=poissons_ratio,
        yield_strength=yield_strength,
        traction_values=traction_values,
    )

    # Return the factor of safety
    return fos
