"""Lorentzian geometric existence certificate for PCRG's four-move history.

The coordinates below embed all four 4-simplices in one Minkowski space with
signature (-,+,+,+).  Every tetrahedron appearing on any intermediate
hypersurface is spacelike, and the outward timelike normals have exactly the
past/future pattern required by the 2->3, 1->4, 2->3, 3->2 sequence.

This certifies a nondegenerate time-oriented Lorentzian geometry.  It does not
choose a complex Regge-action branch or an integration contour.
"""

from itertools import combinations

import numpy as np


MINKOWSKI_METRIC = np.diag((-1.0, 1.0, 1.0, 1.0))

# A reproducible numerical witness.  The time coordinate of a random rationally
# independent seed geometry was enlarged by 5/4 to improve the Lorentzian-volume
# margin while retaining spacelike hypersurfaces.
LORENTZIAN_POINTS = np.array(
    [
        [0.09231901942494321, -0.26623686693716175, -0.9561471332773008, 0.5263766901523332],
        [0.00779087413491134, 0.6641656902648901, 0.6387994228952245, 0.3273332669933098],
        [-0.07177288627508509, 0.00549547741536736, 0.580942220887827, -1.3904473720084403],
        [-0.05108933380926729, -1.1349286507775327, -0.2132660713209228, 0.7650241827240767],
        [0.02584866917855915, 0.3933809084163796, -0.09483793096062454, -0.3057058993858946],
        [-0.00309634265406130, 0.3381234416180574, 0.04450949177579638, 0.07741913152461477],
    ],
    dtype=float,
)

BULK_SIMPLICES = (
    (0, 1, 2, 3, 4),
    (0, 1, 3, 4, 5),
    (1, 2, 3, 4, 5),
    (0, 2, 3, 4, 5),
)

PAST_FACES = (
    {(0, 1, 2, 3), (0, 1, 2, 4)},
    {(0, 1, 3, 4)},
    {(1, 2, 3, 4), (1, 3, 4, 5)},
    {(0, 2, 3, 4), (0, 3, 4, 5), (2, 3, 4, 5)},
)

A7_TO_HISTORY = {0: 3, 1: 4, 2: 0, 3: 1, 4: 2, 5: 5}
A7_FOUR_TWO_SIGN_BITS = np.asarray((0, 0, 1, 1, 1, 1), dtype=int)
INTERNAL_TRIANGLES = ((0, 3, 4), (1, 3, 4), (2, 3, 4), (3, 4, 5))


def minkowski_product(first: np.ndarray, second: np.ndarray) -> float:
    return float(first @ MINKOWSKI_METRIC @ second)


def induced_gram(vertices: tuple[int, ...], points: np.ndarray) -> np.ndarray:
    vectors = np.asarray([points[vertex] - points[vertices[0]] for vertex in vertices[1:]])
    return vectors @ MINKOWSKI_METRIC @ vectors.T


def outward_unit_normal(
    face: tuple[int, int, int, int], opposite: int, points: np.ndarray
) -> np.ndarray:
    origin = points[face[0]]
    tangents = np.asarray([points[vertex] - origin for vertex in face[1:]])
    # If n is a contravariant normal, tangents @ eta @ n = 0.
    normal = np.linalg.svd(tangents @ MINKOWSKI_METRIC)[2][-1]
    if minkowski_product(normal, points[opposite] - origin) > 0.0:
        normal = -normal
    norm_squared = minkowski_product(normal, normal)
    if norm_squared >= 0.0:
        raise ValueError("facet normal is not timelike")
    return normal / np.sqrt(-norm_squared)


def source_mapped_pachner_hessian_entry(
    first_edge: tuple[int, int],
    second_edge: tuple[int, int],
    d_factor: float,
    volumes: np.ndarray,
) -> float:
    """Published factorized Lorentzian 4-2 Hessian in A7 labels."""
    first, second = first_edge
    third, fourth = second_edge
    exponent = (
        A7_FOUR_TWO_SIGN_BITS[first]
        + A7_FOUR_TWO_SIGN_BITS[second]
        + A7_FOUR_TWO_SIGN_BITS[third]
        + A7_FOUR_TWO_SIGN_BITS[fourth]
    )
    sign = (-1.0) ** int(exponent)
    numerator = volumes[first] * volumes[second] * volumes[third] * volumes[fourth]
    return float(-d_factor * sign * numerator / (48.0 * np.prod(volumes)))


def complex_log_minus(value: complex, tolerance: float = 1e-12) -> complex:
    """Principal logarithm with argument -pi on the negative-real branch cut."""
    value = complex(value)
    if value.real < 0.0 and abs(value.imag) < tolerance:
        return complex(np.log(-value.real), -np.pi)
    return complex(np.log(value))


def complex_sqrt_plus(value: complex) -> complex:
    """Principal square root, equal to +i sqrt(r) for a negative real -r."""
    return complex(np.sqrt(complex(value)))


def complex_dihedral_angle(
    simplex: tuple[int, int, int, int, int],
    hinge: tuple[int, int, int],
    points: np.ndarray,
) -> complex:
    """Source-convention theta-plus from outward normals and the log-minus branch.

    The second outward normal is reversed so the vector angle is the interior
    dihedral angle.  On Euclidean data this convention gives minus the usual
    positive interior angle, as required by ``epsilon=2*pi+sum(theta)``.
    """
    opposite = [vertex for vertex in simplex if vertex not in hinge]
    if len(opposite) != 2:
        raise ValueError("hinge must be a triangle in the four-simplex")
    first_face = tuple(vertex for vertex in simplex if vertex != opposite[0])
    second_face = tuple(vertex for vertex in simplex if vertex != opposite[1])
    first = outward_unit_normal(first_face, opposite[0], points)
    second = -outward_unit_normal(second_face, opposite[1], points)
    first_norm = minkowski_product(first, first)
    second_norm = minkowski_product(second, second)
    product = minkowski_product(first, second)
    discriminant = first_norm * second_norm - product**2
    argument = (
        product - 1j * complex_sqrt_plus(discriminant)
    ) / (complex_sqrt_plus(first_norm) * complex_sqrt_plus(second_norm))
    return complex(-1j * complex_log_minus(argument))


def internal_complex_deficits(points: np.ndarray) -> dict:
    angles = {}
    deficits = {}
    lightcone_pi_counts = {}
    for hinge in INTERNAL_TRIANGLES:
        incident_angles = []
        for simplex in BULK_SIMPLICES:
            if set(hinge) <= set(simplex):
                angle = complex_dihedral_angle(simplex, hinge, points)
                angles[(simplex, hinge)] = angle
                incident_angles.append(angle)
        deficits[hinge] = complex(2.0 * np.pi + sum(incident_angles))
        lightcone_pi_counts[hinge] = sum(
            abs(angle.real + np.pi) < 2e-12 for angle in incident_angles
        )
    return {
        "angles": angles,
        "deficits": deficits,
        "lightcone_pi_counts": lightcone_pi_counts,
        "maximum_deficit_residual": max(abs(value) for value in deficits.values()),
    }


def lorentzian_history_certificate(points: np.ndarray = LORENTZIAN_POINTS) -> dict:
    points = np.asarray(points, dtype=float)
    if points.shape != (6, 4):
        raise ValueError("points must have shape (6,4)")

    squared_intervals = {
        edge: minkowski_product(points[edge[1]] - points[edge[0]], points[edge[1]] - points[edge[0]])
        for edge in combinations(range(6), 2)
    }

    simplex_data = []
    face_normals: dict[tuple[int, int, int, int], list[np.ndarray]] = {}
    directed_time_components = []
    tetrahedron_min_eigenvalues = []
    for simplex_index, simplex in enumerate(BULK_SIMPLICES):
        gram = induced_gram(simplex, points)
        eigenvalues = np.linalg.eigvalsh(gram)
        coordinate_matrix = np.asarray(
            [points[vertex] - points[simplex[0]] for vertex in simplex[1:]]
        )
        if not (eigenvalues[0] < 0.0 < eigenvalues[1]):
            raise ValueError("four-simplex is not Lorentzian")

        normals = {}
        for opposite in simplex:
            face = tuple(vertex for vertex in simplex if vertex != opposite)
            face_eigenvalues = np.linalg.eigvalsh(induced_gram(face, points))
            if face_eigenvalues[0] <= 0.0:
                raise ValueError("history tetrahedron is not spacelike")
            tetrahedron_min_eigenvalues.append(float(face_eigenvalues[0]))
            normal = outward_unit_normal(face, opposite, points)
            normals[face] = normal
            face_normals.setdefault(face, []).append(normal)
            desired_sign = -1.0 if face in PAST_FACES[simplex_index] else 1.0
            directed_time = desired_sign * normal[0]
            if directed_time <= 0.0:
                raise ValueError("facet has the wrong past/future time orientation")
            directed_time_components.append(float(directed_time))

        simplex_data.append(
            {
                "simplex": simplex,
                "gram": gram,
                "gram_eigenvalues": eigenvalues,
                "coordinate_determinant": float(np.linalg.det(coordinate_matrix)),
                "normals": normals,
            }
        )

    shared_normal_residuals = {
        face: float(np.linalg.norm(normals[0] + normals[1]))
        for face, normals in face_normals.items()
        if len(normals) == 2
    }
    if len(shared_normal_residuals) != 6:
        raise ValueError("unexpected number of shared tetrahedra")

    homogeneous = np.vstack((np.ones(6), points.T))
    circuit = np.linalg.svd(homogeneous)[2][-1]
    circuit /= np.max(np.abs(circuit))
    if set(np.flatnonzero(circuit > 0.0)) != {0, 1, 2, 5}:
        circuit = -circuit
    if set(np.flatnonzero(circuit > 0.0)) != {0, 1, 2, 5}:
        raise ValueError("affine circuit has the wrong 4-2 Radon partition")

    # Source-mapped A8 Lorentzian 4-2 Hessian.  Relabel the history so that its
    # internal edge (34) is A7 edge (01), then use the published volume/D formula.
    a7_points = np.asarray([points[A7_TO_HISTORY[vertex]] for vertex in range(6)])
    a7_homogeneous = np.vstack((np.ones(6), a7_points.T))
    a7_circuit = np.asarray(
        [
            (-1.0) ** omitted
            * np.linalg.det(np.delete(a7_homogeneous, omitted, axis=1))
            / 24.0
            for omitted in range(6)
        ]
    )
    if a7_circuit[0] < 0.0:
        a7_circuit = -a7_circuit
    if not (np.all(a7_circuit[:2] > 0.0) and np.all(a7_circuit[2:] < 0.0)):
        raise ValueError("relabeled affine circuit does not realize the A7 4-2 partition")
    affine_d_factor = float(
        2.0
        * sum(
            a7_circuit[vertex]
            * minkowski_product(a7_points[vertex], a7_points[vertex])
            for vertex in range(6)
        )
    )
    a7_volumes = []
    for omitted in range(6):
        simplex_points = np.delete(a7_points, omitted, axis=0)
        a7_volumes.append(
            abs(np.linalg.det((simplex_points[1:] - simplex_points[0]).T)) / 24.0
        )
    a7_volumes = np.asarray(a7_volumes)
    # The coordinate circuit convention of A4 is eight times the D used in the
    # primary Lorentzian Hessian formula.  Verify that conversion independently
    # from the source's explicit 4-2 edge-length expression.
    published_d_factor = affine_d_factor / 8.0

    def signed_interval(first: int, second: int) -> float:
        difference = a7_points[first] - a7_points[second]
        return minkowski_product(difference, difference)

    published_d_from_lengths = float(
        -sum(
            (
                signed_interval(0, vertex)
                + signed_interval(1, vertex)
                - signed_interval(0, 1)
            )
            * a7_volumes[vertex]
            for vertex in range(2, 6)
        )
        / 8.0
    )
    radial_hessian = float(
        -published_d_factor
        * a7_volumes[0]
        * a7_volumes[1]
        / (48.0 * np.prod(a7_volumes[2:]))
    )
    input_cross = source_mapped_pachner_hessian_entry(
        (0, 1), (2, 3), published_d_factor, a7_volumes
    )
    output_cross = source_mapped_pachner_hessian_entry(
        (0, 1), (4, 5), published_d_factor, a7_volumes
    )
    pachner_boundary_cross = source_mapped_pachner_hessian_entry(
        (2, 3), (4, 5), published_d_factor, a7_volumes
    )
    separated_transfer = float(-input_cross * output_cross / radial_hessian)
    complex_curvature = internal_complex_deficits(points)

    return {
        "points": points,
        "squared_intervals": squared_intervals,
        "simplex_data": simplex_data,
        "minimum_spacelike_edge_interval": min(squared_intervals.values()),
        "minimum_tetrahedron_gram_eigenvalue": min(tetrahedron_min_eigenvalues),
        "minimum_directed_normal_time_component": min(directed_time_components),
        "shared_normal_residuals": shared_normal_residuals,
        "maximum_shared_normal_residual": max(shared_normal_residuals.values()),
        "affine_circuit": circuit,
        "affine_circuit_residual": float(np.linalg.norm(homogeneous @ circuit)),
        "a7_affine_circuit": a7_circuit,
        "a7_absolute_four_volumes": a7_volumes,
        "lorentzian_affine_d_factor": affine_d_factor,
        "lorentzian_d_factor": published_d_factor,
        "lorentzian_d_from_lengths": published_d_from_lengths,
        "d_normalization_residual": abs(published_d_factor - published_d_from_lengths),
        "source_mapped_radial_hessian": radial_hessian,
        "source_mapped_input_internal_cross": input_cross,
        "source_mapped_internal_output_cross": output_cross,
        "source_mapped_pachner_boundary_cross": pachner_boundary_cross,
        "source_mapped_separated_transfer": separated_transfer,
        "rank_one_transfer_residual": abs(separated_transfer + pachner_boundary_cross),
        "complex_internal_angles": complex_curvature["angles"],
        "complex_internal_deficits": complex_curvature["deficits"],
        "lightcone_pi_counts": complex_curvature["lightcone_pi_counts"],
        "maximum_complex_deficit_residual": complex_curvature[
            "maximum_deficit_residual"
        ],
        "gaussian_steepest_descent_phase": float(np.sign(radial_hessian) * np.pi / 4.0),
    }


if __name__ == "__main__":
    certificate = lorentzian_history_certificate()
    print(
        "minimum_spacelike_edge_interval="
        f"{certificate['minimum_spacelike_edge_interval']:.16g}"
    )
    print(
        "minimum_tetrahedron_gram_eigenvalue="
        f"{certificate['minimum_tetrahedron_gram_eigenvalue']:.16g}"
    )
    for entry in certificate["simplex_data"]:
        print(
            f"simplex={entry['simplex']} determinant={entry['coordinate_determinant']:.16g} "
            f"gram_eigenvalues={entry['gram_eigenvalues']}"
        )
    print(
        "minimum_directed_normal_time_component="
        f"{certificate['minimum_directed_normal_time_component']:.16g}"
    )
    print(
        "maximum_shared_normal_residual="
        f"{certificate['maximum_shared_normal_residual']:.3e}"
    )
    print(f"affine_circuit={certificate['affine_circuit']}")
    print(f"affine_circuit_residual={certificate['affine_circuit_residual']:.3e}")
    print(
        "lorentzian_affine_D="
        f"{certificate['lorentzian_affine_d_factor']:.16g}"
    )
    print(f"lorentzian_source_D={certificate['lorentzian_d_factor']:.16g}")
    print(f"D_normalization_residual={certificate['d_normalization_residual']:.3e}")
    print(
        "source_mapped_radial_hessian="
        f"{certificate['source_mapped_radial_hessian']:.16g}"
    )
    print(
        "gaussian_steepest_descent_phase="
        f"{certificate['gaussian_steepest_descent_phase']:.16g}"
    )
    print(
        "input_internal_cross="
        f"{certificate['source_mapped_input_internal_cross']:.16g}"
    )
    print(
        "internal_output_cross="
        f"{certificate['source_mapped_internal_output_cross']:.16g}"
    )
    print(
        "lorentzian_separated_transfer="
        f"{certificate['source_mapped_separated_transfer']:.16g}"
    )
    print(f"rank_one_transfer_residual={certificate['rank_one_transfer_residual']:.3e}")
    print(
        "maximum_complex_deficit_residual="
        f"{certificate['maximum_complex_deficit_residual']:.3e}"
    )
    print(f"lightcone_pi_counts={certificate['lightcone_pi_counts']}")
