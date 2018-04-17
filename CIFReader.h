//
// Created by Jon on 21/08/2015.
//

#ifndef XYZ_CIFREADER_H
#define XYZ_CIFREADER_H

// for file input
#include <string>
#include <fstream>
#include <streambuf>
#include <regex>
#include <sstream>
#include <vector>

#include "Utilities.h"
#include "AtomSite.h"
#include "SymmetryOperation.h"
#include "UnitCell.h"
#include "CellGeometry.h"

class CIFReader
{
public:
    // constructor that takes filename
    CIFReader(std::string filePath);

    // method to return a class instance that will contain the unit cell information
    UnitCell getUnitCell() { return UnitCell(cell, atomsites); }

private:
    // class instance to hold unit cell information
    std::vector<Symmetry> symmetrylist;

    std::vector<AtomSite> atomsites;

    CellGeometry cell;

    // methods to read
    // 1. atom positions
    // 2. symmetry
    // 3. basis geometry
    void readAtomPositions(const std::string& input);

    void readSymmetryOperations(const std::string& input);

    void readCellGeometry(const std::string& input);

};


#endif //XYZ_CIFREADER_H
