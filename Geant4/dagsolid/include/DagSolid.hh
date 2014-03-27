//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration and of QinetiQ Ltd,   *
// * subject to DEFCON 705 IPR conditions.                            *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: DagSolid.hh,v 1.10 2010/12/10 16:30:13 gunter Exp $
// GEANT4 tag $Name: geant4-09-05 $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:              DagSolid.hh
//
// Date:                20/12/2010
// Author:              M. C. Han, C. H. Kim, J. H. Jeong, Y. S. Yeom, S. Kim, 
//                      Paul. P. H. Wilson, J. Apostolakis
// Organisation:        Hanyang Univ., KR
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 31 October 2010, J. H. Jeong, Hanyang Univ., KR
//  - Created.
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



///////////////////////////////////////////////////////////////////////////////
#ifndef DagSolid_hh
#define DagSolid_hh 1

#include "G4TessellatedSolid.hh"
#include "G4VGraphicsScene.hh"
#include "G4VPVParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "globals.hh"

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "DagMC.hpp"
using namespace moab;


class DagSolid : public G4TessellatedSolid
{
  public:  // with description

    DagSolid ();
    DagSolid (const G4String &name, DagMC* dagmc, int volID);
    virtual ~DagSolid ();


    // Mandatory Functions
    
    virtual EInside Inside (const G4ThreeVector &p) const;

    virtual G4ThreeVector SurfaceNormal (const G4ThreeVector &p) const;
    virtual G4double DistanceToIn(const G4ThreeVector &p,
                                  const G4ThreeVector &v) const;
    virtual G4double DistanceToIn(const G4ThreeVector &p) const;
    virtual G4double DistanceToOut(const G4ThreeVector &p,
                                   const G4ThreeVector &v,
                                   const G4bool calcNorm=false,
                                   G4bool *validNorm=0, G4ThreeVector *n=0) const;
    virtual G4double DistanceToOut (const G4ThreeVector &p) const;
    virtual G4GeometryType GetEntityType () const;

    virtual G4bool CalculateExtent(const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, G4double& pMax) const;
    virtual std::ostream &StreamInfo(std::ostream &os) const;



//    virtual void ComputeDimensions (G4VPVParameterisation* p, const G4int n,
//                                  const G4VPhysicalVolume* pRep);
 

    virtual G4double GetCubicVolume ();
    virtual G4double GetSurfaceArea ();
    G4double      GetMinXExtent () const;
    G4double      GetMaxXExtent () const;
    G4double      GetMinYExtent () const;
    G4double      GetMaxYExtent () const;
    G4double      GetMinZExtent () const;
    G4double      GetMaxZExtent () const;
  // Functions for visualization
 
    virtual void  DescribeYourselfTo (G4VGraphicsScene& scene) const;



  public:  // without description

    DagSolid(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.



  protected:  // with description
 
    void DeleteObjects ();
    void CopyObjects (const DagSolid &s);


  private:


    G4GeometryType           geometryType;
    G4double                 cubicVolume;
    G4double                 surfaceArea;
    G4double                 xMinExtent;
    G4double                 xMaxExtent;
    G4double                 yMinExtent;
    G4double                 yMaxExtent;
    G4double                 zMinExtent;
    G4double                 zMaxExtent;
  

    G4String Myname;
     
    DagMC* fdagmc;
    G4int fvolID;
    EntityHandle fvolEntity;

    mutable EntityHandle Last_sulf_hit;
    mutable G4int nVertices;
    EntityHandle My_sulf_hit;    
    
};

#endif
