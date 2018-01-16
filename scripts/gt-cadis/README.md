GT-CADIS
========

This workflow was initially developed as part of the PhD research of Elliott
Biondo. The theory and background can be found in the dissertation [1] produced as
part of that work.

Development plan
----------------

The current development plan is outlined below.  Much of the capability
described here is available in some form or approximation in the archive [2] made
available by Biondo on Figshare.

The following steps describe a complete GT-CADIS workflow.  It is designed to
be highly modular with the potential to combine pieces into larger scripts as
users and developers find convenient.

1. Discretize geometry for deterministic transport
1. Prepare for adjoint photon transport

    1. combine discretized geom with adjoint source definition
    1. write PARTISN input file

1. Run PARTISN

1. Post-process adjoint photon transport
    1. read ATFLUX adj photon flux file
    1. write adjoint photon flux in convenient format (H5M, VTK)

1. Prepare ALARA files to generate T matrices, one per material, per irradiation history (option for future feature: different T matrices with various non-SNILB corrections)
    1. requires a flux spectrum of some kind - could be multiple (ideally small N)
        1. reference spectrum
        1. (option) deterministic spectrum from this geometry (which voxel??)
    1. prepare ALARA input file
        1. for each material
        1. for each input spectrum
        1. for each neutron group
        1. one "volume" with a one-group spectrum  (robustness: what if one group has 0 flux?)
        1. also one for the total neutron flux
    1. write ALARA input files

1. Run ALARA

1. Calculate T matrices and perform SNILB tests
    1. process ALARA output file
        1. test SNILB with multiple spectrums:  S(spectrum) =? sum(one-group sources)
        1. calculate T from one of the N fluxes used to test SNILB
    1. write SNILB eta values
    1. write T matrices to disk (option)

1. Calculate adjoint neutron source
    1. read adjoint photon flux from ATFLUX or H5M/VTK
    1. read T matrices
    1. calculate adjoint neutron source
    1. write adjoint neutron source

1. Prepare adjoint neutron transport
    1. combine discretized geom with adjoint source definition
    1. write PARTISN input file

1. Run PARTISN

1. Post-process adjoint neutron transport
    1. read ATFLUX adj neutron flux file
    1. write adjoint neutron flux in convenient format (H5M, VTK)

1. Neutron CADIS
    1. read adjoint neutron fluxes from ATFLUX or H5M/VTK
    1. read unbiased forward neutron source
    1. calculate source biases and WW lower bounds
    1. write biased forward neutron source
    1. write WWINP and H5M/VTK format

1. Perform neutron transport with DAG-MCNP

1. R2S-CADIS

File inventory:

Inputs
---------
* DAGMC geometry w/ detector
* Irradiation history w/ cooling times (R2S step 2 information)
* unbiased forward neutron source mesh
* DAG-MCNP input files

Generated
---------------
* T matrices - one per material (format:???) [viz]
* SNILB eta results - N per material - output to screen? can be saved by user?
* Partisn input file for adjoint photon transport
* Adj photon flux files (ATFLUX, H5M/VTK) [viz]
* Partisn input file for adjoint neutron transport
* Adj neutron flux files (ATFLUX, H5M/VTK) [viz]
* ALARA input files
* ALARA output files
* biased forward neutron source mesh [viz]
* weight window lower bounds (WWINP, H5M/VTK) [viz]


[1]: http://depot.library.wisc.edu/repository/fedora/1711.dl:ITANHEGGRPM338Z/datastreams/REF/content

[2]: https://figshare.com/articles/Supporting_files_for_Transmutation_Approximations_for_the_Application_of_Hybrid_Monte_Carlo_Deterministic_Neutron_Transport_to_Shutdown_Dose_Rate_Analysis_/3546432
