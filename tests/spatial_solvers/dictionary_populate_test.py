"""This file creates a sample dictionary to test the spatial_solver code's."""
import numpy as np

cell_widths = np.array([0.25, 0.25, 0.25, 0.25], "float64", order="F")


def populate_simple(solverin, solvertypein):
    solver_dict = {
        "solver": solverin,
        "solver_type": solvertypein,
        "spatial_order": 1,
        "angular_quadrature_order": 4,
        "angular_quadrature_type": 1,
        "nodes_xyz": [4, 4, 4],
        "num_groups": 1,
        "num_materials": 1,
        "x_cells_widths": cell_widths,
        "y_cells_widths": cell_widths,
        "z_cells_widths": cell_widths,
        #'x_cells_widths':[0.25, 0.25, 0.25, 0.25],
        #'y_cells_widths':[0.25, 0.25, 0.25, 0.25],
        #'z_cells_widths':[0.25, 0.25, 0.25, 0.25],
        "x_boundry_conditions": [2, 2],
        "y_boundry_conditions": [2, 2],
        "z_boundry_conditions": [2, 2],
        #'material_id': np.array([[[1]*4]*4]*4, 'float64'),
        "material_id": [[[1] * 4] * 4] * 4,
        "quadrature_file": "quad_file",
        "xs_file": "spatial_solvers/xs",
        "source_input_file": "spatial_solvers/src_4.dat",
        #'source_input_file':'src_non_binary',
        "bc_input_file": "spatial_solvers/bc_4.dat",
        "flux_output_file": "spatial_solvers/phi_4.ahot",
        "convergence_criterion": 1.0e-12,
        "max_iterations": 6000,
        "moments_converged": 0,
        "converge_tolerence": 1.0e-10,
    }

    if solverin == "SCTSTEP":
        solver_dict["bc_input_file"] = "spatial_solvers/bc_4_sct.dat"
    return solver_dict


def populate_simple_with_warnings(solverin):
    solver_dict = {
        "solver": solverin,
        #'solver_type':solvertypein,
        #'spatial_order':1,
        #'angular_quadrature_order':4,
        #'angular_quadrature_type':1,
        "nodes_xyz": [4, 4, 4],
        "num_groups": 1,
        "num_materials": 1,
        "x_cells_widths": cell_widths,
        "y_cells_widths": cell_widths,
        "z_cells_widths": cell_widths,
        "x_boundry_conditions": [2, 2],
        "y_boundry_conditions": [2, 2],
        "z_boundry_conditions": [2, 2],
        "material_id": [[[1] * 4] * 4] * 4,
        "quadrature_file": "spatial_solvers/quad_file",
        "xs_file": "spatial_solvers/xs",
        "source_input_file": "spatial_solvers/src_4.dat",
        "bc_input_file": "spatial_solvers/bc_4.dat",
        "flux_output_file": "spatial_solvers/phi_4.ahot",
        #'convergence_criterion':1.e-12,
        #'max_iterations':6000,
        #'moments_converged':0,
        #'converge_tolerence':1.e-10
    }
    return solver_dict


def populate_intermediate_1(solverin, solvertypein):
    solver_dict = {
        "solver": solverin,
        "solver_type": solvertypein,
        "spatial_order": 1,
        "angular_quadrature_order": 4,
        "angular_quadrature_type": 1,
        "nodes_xyz": [4, 4, 4],
        "num_groups": 1,
        "num_materials": 2,
        "x_cells_widths": cell_widths,
        "y_cells_widths": cell_widths,
        "z_cells_widths": cell_widths,
        "x_boundry_conditions": [2, 2],
        "y_boundry_conditions": [2, 2],
        "z_boundry_conditions": [2, 2],
        "material_id": [
            [
                [1, 2, 1, 2],
                [2, 1, 2, 1],
                [1, 2, 1, 2],
                [2, 1, 2, 1],
            ],
            [
                [2, 1, 2, 1],
                [1, 2, 1, 2],
                [2, 1, 2, 1],
                [1, 2, 1, 2],
            ],
            [
                [1, 2, 1, 2],
                [2, 1, 2, 1],
                [1, 2, 1, 2],
                [2, 1, 2, 1],
            ],
            [
                [2, 1, 2, 1],
                [1, 2, 1, 2],
                [2, 1, 2, 1],
                [1, 2, 1, 2],
            ],
        ],
        "quadrature_file": "quad_file",
        "xs_file": "spatial_solvers/xs_alternating",
        "source_input_file": "spatial_solvers/src_4.dat",
        "bc_input_file": "spatial_solvers/bc_4.dat",
        "flux_output_file": "spatial_solvers/phi_4.ahot",
        "convergence_criterion": 1.0e-12,
        "max_iterations": 6000,
        "moments_converged": 0,
        "converge_tolerence": 1.0e-10,
    }
    if solverin == "SCTSTEP":
        solver_dict["bc_input_file"] = "spatial_solvers/bc_4_sct.dat"
    return solver_dict
