"""This module accesses the AHOTN, DGFEM and SCT-STEP flavored nuetron transport solvers.
"""
from __future__ import division

import sys
import os
import time
import numpy as np

# Solver imports
import pyne.transport_spatial_methods as transport_spatial_methods


def solve(inputdict_unchecked):
    """
    Input Dictionary Required Key Pair Values
    +-----------------------+--------------+----------------------------------+--------------------------+
    |          Key          |     Type     |           Description            |         Example          |
    +-----------------------+--------------+----------------------------------+--------------------------+
    |         solver        |    String    |              Solver              |          AHOTN           |
    |      solver_type      |    String    |           Solver type            |            LN            |
    |     spatial_order     |     Int      |     Spatial expansion order      |            1             |
    |    quadrature_order   |     Int      |         Angles per Octet         |            4             |
    |    quadrature_type    |     Int      |        Quadrature Scheme         |            1             |
    |       nodes_xyz       |  Int Array   | Spatial nodes in XY&Z Directions |        [4, 4, 4]         |
    |       num_groups      |     Int      |     Number of Energy Groups      |            1             |
    |     num_materials     |     Int      |       Number of Materials        |            1             |
    |     x_cells_widths    | Double Array |           Cell X Size            | [0.25, 0.25, 0.25, 0.25] |
    |     y_cells_widths    | Double Array |           Cell Y Size            | [0.25, 0.25, 0.25, 0.25] |
    |     z_cells_widths    | Double Array |           Cell Z Size            | [0.25, 0.25, 0.25, 0.25] |
    |  x_boundry_conditions |  Int Array   |       X Boundry Conditions       |          [2,2]           |
    |  y_boundry_conditions |  Int Array   |       Y Boundry Conditions       |          [2,2]           |
    |  z_boundry_conditions |  Int Array   |       Z Boundry Conditions       |          [2,2]           |
    |      material_id      | 3d Int Array |     Material Mesh Type Info      |           n/a            |
    |       quad_file       |    String    |   Only needed for quad type 2    |        quad_file         |
    |        xs_file        |    String    |     Cross Section Data File      |         xs_file          |
    |   source_input_file   |    String    |      Source Input Info File      |         src.dat          |
    |     bc_input_file     |    String    |   Boundry Condition Input File   |          bc.dat          |
    |    flux_output_file   |    String    |       File for Flux Write        |         flux.out         |
    | convergence_criterion |    float     |      Convergence Criterion       |          1.e-5           |
    |   converge_tolerance  |    float     |            Tolerance             |          1.e-10          |
    |     max_iterations    |     Int      |     Max Iterations for Sweep     |          10000           |
    |   moments_converged   |     Int      | Moments Converged Upon Per Quad  |            0             |
    +-----------------------+--------------+----------------------------------+--------------------------+

    Output Dictionary Key Pair Values
    +------------+--------------+-----------------------------+
    |    Key     |     Type     |         Description         |
    +------------+--------------+-----------------------------+
    |   solver   |    String    |            Solver           |
    |    flux    | Double Array |  Flux Solution Output Array |
    |  success   |     Int      | Solver Success Code (1 yes) |
    | time_start |    Double    |      System Start Time      |
    | total_time |    Double    |    Total Solver Run Time    |
    | print_time |    Double    | Time Taken to Print Results |
    | error_msg  |    String    |  Error Message (if failure) |
    +------------+--------------+-----------------------------+'

    For more detailed information about the input & output key-pair values, see:
        http://pyne.io/pyapi/spatialsolver.html

    """
    inputdict = _dict_complete(inputdict_unchecked)
    fortran_returns = transport_spatial_methods.main(
        "test title in",
        inputdict["solver"],
        inputdict["solver_type"],
        inputdict["spatial_order"],
        # inputdict['spatial_method'],
        inputdict["angular_quadrature_order"],
        inputdict["angular_quadrature_type"],
        inputdict["nodes_xyz"][0],
        inputdict["nodes_xyz"][1],
        inputdict["nodes_xyz"][2],
        inputdict["num_groups"],
        inputdict["num_materials"],
        inputdict["x_cells_widths"],
        inputdict["y_cells_widths"],
        inputdict["z_cells_widths"],
        #       np.asfortranarray(inputdict['y_cells_widths']),
        #       np.asfortranarray(inputdict['z_cells_widths']),
        inputdict["x_boundry_conditions"][0],
        inputdict["x_boundry_conditions"][1],
        inputdict["y_boundry_conditions"][0],
        inputdict["y_boundry_conditions"][1],
        inputdict["z_boundry_conditions"][0],
        inputdict["z_boundry_conditions"][1],
        inputdict["material_id"],
        inputdict["quadrature_file"],
        inputdict["xs_file"],
        inputdict["source_input_file"],
        inputdict["bc_input_file"],
        inputdict["flux_output_file"],
        inputdict["convergence_criterion"],
        inputdict["max_iterations"],
        inputdict["moments_converged"],
        inputdict["converge_tolerence"],
        inputdict["max_mom_printed"],
        inputdict["moment_sum_flag"],
        inputdict["mom_at_a_pt_flag"],
        inputdict["quad_flux_print_flag"],
    )

    solver_output = inputdict
    solver_output["flux"] = fortran_returns[6].tolist()
    error_code = fortran_returns[7]
    tsolve = fortran_returns[8]
    ttosolve = fortran_returns[9]
    tend = fortran_returns[10]

    if error_code == 0:
        solver_output["success"] = 1
        solver_output["time_start"] = tsolve
        solver_output["total_time"] = tsolve - ttosolve
        solver_output["print_time"] = tend - ttosolve
        solver_output["error_msg"] = 0
    else:
        solver_output["success"] = 0
        solver_output["error_msg"] = _error_to_string(error_code)
        print(solver_output["error_msg"])
    return solver_output


def _error_to_string(error_code):
    err_dictionary = {
        1001: "ERROR: Lambda must be equal to one.",
        1002: "ERROR: Illegal value for qdord. Must be greater than zero.",
        1003: "ERROR: Illegal value for the quadrature order. Even #s only.",
        1004: "ERROR: illegal entry for the qdfile name.",
        1005: "ERROR: illegal entry for lambda. Must be zero or greater.",
        1007: "ERROR: Illegal number of x cells. Must be positive.",
        1008: "ERROR: Illegal number of y cells. Must be positive.",
        1009: "ERROR: Illegal number of z cells. Must be positive.",
        1010: "ERROR: Illegal number of energy groups. Must be positive.",
        1011: "ERROR: Illegal number of materials. Must be positive.",
        1012: "ERROR: Illegal x cell dimension, dx. Must be positive.",
        1013: "ERROR: Illegal y cell dimension, dy. Must be positive.",
        1014: "ERROR: Illegal z cell dimension, dz. Must be positive.",
        1015: "ERROR: Illegal value in material map. Must be in [1, #materials].",
        1016: "ERROR: Illegal lower x BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1017: "ERROR: Illegal lower y BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1018: "ERROR: Illegal lower z BC. Must be 0-Vac, 1-Refl or 2-Fixed.",
        1019: "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed.",
        1020: "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed.",
        1021: "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed.",
        1022: "ERROR: Illegal upper x BC. Must be 0-Vac or 2-Fixed.",
        1023: "ERROR: Illegal upper y BC. Must be 0-Vac or 2-Fixed.",
        1024: "ERROR: Illegal upper z BC. Must be 0-Vac or 2-Fixed.",
        1025: "ERROR: Illegal convergence criterion. Must be positive.",
        1026: "ERROR: Illegal max inner iterations, itmx. Must be positive.",
        1027: "ERROR: Illegal tolerance setting, tolr. Must be positive.",
        1028: "ERROR: Illegal value for moments to converge, iall. Must be in [0, lambda].",
        1029: "ERROR: Illegal value for solution check flag, iall. iall is 0 (skip check) or positive integer",
        1030: "ERROR: Illegal value for solution check tolerance, tchk. Must be positive.",
        1031: "ERROR: Illegal value for direction cosine. Must be entered positive, less than 1.",
        1032: "ERROR: Illegal value for weight. Must be entered positive, less than 0.125.",
        1033: "ERROR: a discrete ordinate has length not in range 0.99<1.00<1.01 based on mu, eta, xi values.",
        1034: "ERROR: Illegal value for max moment to print, momp. Must be between 0 and lambda.",
        1035: "ERROR: Illegal value for flag for moment summing. Must be 0 for off or 1 for on.",
        1036: "ERROR: Illegal value for flag for moment sum at non-center point. Must be 0/1=off/on.",
        1037: "ERROR: Illegal value for flag for printing average flux of the four quadrants. Must be 0/1 = off/on.",
    }
    return err_dictionary[error_code]


def _dict_complete(inputdict):

    formatted_dict = {}
    try:
        if (
            (inputdict["solver"] == "AHOTN")
            or (inputdict["solver"] == "DGFEM")
            or (inputdict["solver"] == "SCTSTEP")
        ):
            formatted_dict["solver"] = inputdict["solver"]
        else:
            assert 0 == 1, "solver does not exist"
    except:
        assert 0 == 1, "solver key does not exist"
    try:
        formatted_dict["solver_type"] = inputdict["solver_type"]
    except:
        if inputdict["solver"] == "AHOTN":
            formatted_dict["solver_type"] = "LN"
        elif inputdict["solver"] == "DGFEM":
            formatted_dict["solver_type"] = "LD"
    try:
        formatted_dict["spatial_order"] = inputdict["spatial_order"]
    except:
        formatted_dict["spatial_order"] = 1
        warn(warning_msg + " spatial_order value of 1")
    try:
        formatted_dict["angular_quadrature_order"] = inputdict[
            "angular_quadrature_order"
        ]
    except:
        formatted_dict["angular_quadrature_order"] = 4
        warn(warning_msg + " angular_quadrature_order value of 4")
    try:
        formatted_dict["angular_quadrature_type"] = inputdict["angular_quadrature_type"]
    except:
        formatted_dict["qangular_uadrature_type"] = 1
        warn(warning_msg + " angular_quadrature_type value of 1")
    assert "nodes_xyz" in inputdict, "nodes_xyz key not in dict"
    formatted_dict["nodes_xyz"] = inputdict["nodes_xyz"]
    assert "num_groups" in inputdict, "num_groups key not in dict"
    formatted_dict["num_groups"] = inputdict["num_groups"]
    assert "num_materials" in inputdict, "num_materials key not in dict"
    formatted_dict["num_materials"] = inputdict["num_materials"]
    assert "x_cells_widths" in inputdict, "x_cells_widths not in dict"
    formatted_dict["x_cells_widths"] = inputdict["x_cells_widths"]
    assert "y_cells_widths" in inputdict, "y_cells_widths not in dict"
    formatted_dict["y_cells_widths"] = inputdict["y_cells_widths"]
    assert "z_cells_widths" in inputdict, "z_cells_widths not in dict"
    formatted_dict["z_cells_widths"] = inputdict["z_cells_widths"]
    assert "x_boundry_conditions" in inputdict, "x_boundry_conditions not in dict"
    formatted_dict["x_boundry_conditions"] = inputdict["x_boundry_conditions"]
    assert "y_boundry_conditions" in inputdict, "y_boundry_conditions not in dict"
    formatted_dict["y_boundry_conditions"] = inputdict["y_boundry_conditions"]
    assert "z_boundry_conditions" in inputdict, "z_boundry_conditions not in dict"
    formatted_dict["z_boundry_conditions"] = inputdict["z_boundry_conditions"]
    assert "material_id" in inputdict, "material_id not in dict"
    formatted_dict["material_id"] = inputdict["material_id"]
    assert "quadrature_file" in inputdict, "quadrature_file not in dict"
    formatted_dict["quadrature_file"] = inputdict["quadrature_file"]
    assert "xs_file" in inputdict, "xs_file not in dict"
    formatted_dict["xs_file"] = inputdict["xs_file"]
    assert "source_input_file" in inputdict, "source_input_file not in dict"
    formatted_dict["source_input_file"] = inputdict["source_input_file"]
    assert "bc_input_file" in inputdict, "bc_input_file not in dict"
    formatted_dict["bc_input_file"] = inputdict["bc_input_file"]
    assert "flux_output_file" in inputdict, "flux_output_file not in dict"
    formatted_dict["flux_output_file"] = inputdict["flux_output_file"]
    formatted_dict["max_mom_printed"] = 0
    formatted_dict["moment_sum_flag"] = 0
    formatted_dict["mom_at_a_pt_flag"] = 0
    formatted_dict["quad_flux_print_flag"] = 0
    try:
        formatted_dict["convergence_criterion"] = inputdict["convergence_criterion"]
    except:
        formatted_dict["convergence_criterion"] = 1e-12
        warn(warning_msg + " convergence_criterion value of 1e-12")
    try:
        formatted_dict["max_iterations"] = inputdict["max_iterations"]
    except:
        formatted_dict["max_iterations"] = 6000
        warn(warning_msg + " max_iterations value of 6000")
    try:
        formatted_dict["moments_converged"] = inputdict["moments_converged"]
    except:
        formatted_dict["moments_converged"] = 0
        warn(warning_msg + " moments_converged value of 0")
    try:
        formatted_dict["converge_tolerence"] = inputdict["converge_tolerence"]
    except:
        formatted_dict["converge_tolerence"] = 1e-10
        warn(warning_msg + " converge_tolerence value of 1e-10")
    return formatted_dict
