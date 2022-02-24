"""This file creates a sample dictionary to test the spatial_solver code's."""


def populate(solverin, solvertypein):
    solver_dict = {
        "solver": solverin,
        "solver_type": solvertypein,
        "spatial_order": 1,
        "spatial_method": 0,
        "angular_quadrature_order": 4,
        "angular_quadrature_type": 1,
        "nodes_xyz": [4, 4, 4],
        "num_groups": 1,
        "num_materials": 1,
        "x_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "y_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "z_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "x_boundry_conditions": [2, 2],
        "y_boundry_conditions": [2, 2],
        "z_boundry_conditions": [2, 2],
        "material_id": [[[1] * 4] * 4] * 4,
        "quadrature_file": "quad_file",
        "xs_file": "xs",
        "source_input_file": "src_4.dat",
        "bc_input_file": "bc_4.dat",
        "flux_output_file": "phi_4.ahot",
        "convergence_criterion": 1.0e-12,
        "max_iterations": 6000,
        "moments_converged": 0,
        "converge_tolerence": 1.0e-10,
        "ichk": 0,
        "ichk_tolerence": 1.0e-14,
    }
    if solverin == "SCTSTEP":
        solver_dict["bc_input_file"] = "bc_4_sct.dat"
    return solver_dict


def populate_with_warnings(solverin):
    solver_dict = {
        "solver": solverin,
        #'solver_type':solvertypein,
        #'spatial_order':1,
        #'spatial_method':0,
        #'angular_quadrature_order':4,
        #'angular_quadrature_type':1,
        "nodes_xyz": [4, 4, 4],
        "num_groups": 1,
        "num_materials": 1,
        "x_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "y_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "z_cells_widths": [0.25, 0.25, 0.25, 0.25],
        "x_boundry_conditions": [2, 2],
        "y_boundry_conditions": [2, 2],
        "z_boundry_conditions": [2, 2],
        "material_id": [[[1] * 4] * 4] * 4,
        "quadrature_file": "quad_file",
        "xs_file": "xs",
        "source_input_file": "src_4.dat",
        "bc_input_file": "bc_4.dat",
        "flux_output_file": "phi_4.ahot",
        #'convergence_criterion':1.e-12,
        #'max_iterations':6000,
        #'moments_converged':0,
        #'converge_tolerence':1.e-10,
        #'ichk':0,
        #'ichk_tolerence':1.e-14
    }
    return solver_dict
