//
// Created by Jon on 21/08/2015.
//

#ifndef XYZ_UNITCELL_H
#define XYZ_UNITCELL_H

#include "CellGeometry.h"
#include "AtomSite.h"


class UnitCell {
public:
    // constructor that accepts lengths and angles of unit cell
    UnitCell(CellGeometry geo, std::vector<AtomSite> sites) : geometry(geo), atoms(sites) {}

    CellGeometry getCellGeometry(){return geometry;}
    std::vector<AtomSite> getAtoms(){return atoms;}

private:
    CellGeometry geometry;

    std::vector<AtomSite> atoms;
//    CellGeometry geometry;
//
//    std::vector<AtomSite> atoms;

    // method to convert lengths and angles to 3 vectors
};


#endif //XYZ_UNITCELL_H
