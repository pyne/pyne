DagSolid Geant4 Geometry Class
==========================================================

The Direct Acclerated Geometry Monte Carlo (DAGMC) Toolkit is an
interface (`DAGMC Source
<http://trac.mcs.anl.gov/projects/ITAPS/browser/MOAB/trunk/tools/dagmc>`_)
to the `MOAB mesh database
<http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_ that provides the
methods necessary for ray tracing on a CAD-based geometric model.

This repository provides a DAG implimentation for Geant4 (`Geant4 source
<http://geant4.cern.ch/>`_). The Geant4 (G4) implementation is different to other Monte Carlo (MC)
codes in that is specifically a toolkit and thus dictates how the DAG
components are integrated. 

Geometry Interface
-------------------

The DagSolid DAG interface is disributed as a source code which produces a library, it was determined
that this was the most "Geant-like" way of interfacing with G4. The fundamental geometry subroutines
are found in `src/DagSolid.cc` and their associated prototypes are found in `include/DagSolid.hh`. When 
appropriately included in a G4 project (by using the DagSolid methods and includes) users can 
include new DAG volumes in a very G4-like way, for example

1) We must instanciate a new DAG instance
 DagMC* dagmc = DagMC::instance(); // create dag instance
2) Load the h5m file
  const char* h5mfilename = "test_geom.h5m";
3) Having succesfully loaded the file, create a new DagSolid component
  DagSolid* vol_1 = new DagSolid("vol_1",dagmc,1);
4) Thus now DagSolid volumes behave as any other volume, placing logical and physical volumes

  G4LogicalVolume* vol_1_log = new G4LogicalVolume(vol_1,Vacuum,"vol_1_log",0,0,0);
  G4PVPlacement* vol_1_phys = new G4PVPlacement(0,G4ThreeVector(0*cm,0*cm,0*cm),vol_1_log,
                                   "vol_1_phys",world_volume_log,false,0);
5) Interact with these volumes as you would any other G4 volume.


Automated Geant4 Builder
----------------------------------

For those less familiar with G4 and/or those who have a complete DAG model to run, there is an
automated G4 builder which will create a full hierarchy of files including transposing of the 
University of Wisconsin Unified Workflow (UWUW) components already tagged into the h5m geometry.
Given the highly adaptable nature of G4 it is impossible to accomodate every users requirements,
and thus it is recommended to use the G4 builder script as a helper.

