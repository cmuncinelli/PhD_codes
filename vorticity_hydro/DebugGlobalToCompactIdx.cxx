#include <iostream>
#include <stdexcept>
#include <cassert>
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

// ==== Your mapping functions ====

// Forward: ROOT global bin -> compact [0..N-1]
inline int mapGlobalToCompact(const TH1& h, int globalBin){
    int ix, iy, iz;
    h.GetBinXYZ(globalBin, ix, iy, iz);

    int nx = h.GetNbinsX();
    int ny = h.GetNbinsY();
    int nz = h.GetNbinsZ();

    if (ix < 1 || ix > nx)
        throw std::out_of_range("mapGlobalToCompact: X bin under/overflow");
    if (h.GetDimension() >= 2 && (iy < 1 || iy > ny))
        throw std::out_of_range("mapGlobalToCompact: Y bin under/overflow");
    if (h.GetDimension() == 3 && (iz < 1 || iz > nz))
        throw std::out_of_range("mapGlobalToCompact: Z bin under/overflow");

    int cx = ix - 1;
    int cy = (h.GetDimension() >= 2) ? iy - 1 : 0;
    int cz = (h.GetDimension() == 3) ? iz - 1 : 0;

    return cx + nx * (cy + ny * cz);
}

// Inverse: compact [0..N-1] -> ROOT global bin
inline int mapCompactToGlobal(const TH1& h, int compactIndex){
    int nx = h.GetNbinsX();
    int ny = h.GetNbinsY();
    int nz = h.GetNbinsZ();

    int cx = 0, cy = 0, cz = 0;

    if (h.GetDimension() == 3){
        cz = compactIndex / (nx * ny);
        int rem = compactIndex % (nx * ny);
        cy = rem / nx;
        cx = rem % nx;
    }
    else if (h.GetDimension() == 2){
        cy = compactIndex / nx;
        cx = compactIndex % nx;
    }
    else if (h.GetDimension() == 1){
        cx = compactIndex;
    }
    else{
        throw std::logic_error("mapCompactToGlobal: unsupported dimension");
    }

    int ix = cx + 1;
    int iy = (h.GetDimension() >= 2) ? cy + 1 : 0;
    int iz = (h.GetDimension() == 3) ? cz + 1 : 0;

    return h.GetBin(ix, iy, iz);
}

// ==== Test harness ====
template<typename H>
void testHistogram(H& h, const char* name){
    int nBinsTotal = h.GetNbinsX() * std::max(1,h.GetNbinsY()) * std::max(1,h.GetNbinsZ());

    std::cout << "Testing " << name << " with " << nBinsTotal << " bins..." << std::endl;

    // Loop over all global bins
    for (int globalBin = 1; globalBin <= h.GetNcells(); ++globalBin){
        int ix, iy, iz;
        h.GetBinXYZ(globalBin, ix, iy, iz);

        // Skip under/overflow bins
        if (ix < 1 || ix > h.GetNbinsX()) continue;
        if (h.GetDimension() >= 2 && (iy < 1 || iy > h.GetNbinsY())) continue;
        if (h.GetDimension() == 3 && (iz < 1 || iz > h.GetNbinsZ())) continue;

        // Map global -> compact -> global
        int compact = mapGlobalToCompact(h, globalBin);
        int globalBack = mapCompactToGlobal(h, compact);

        if (globalBin != globalBack){
            std::cerr << "Mismatch in " << name
                      << ": globalBin=" << globalBin
                      << " -> compact=" << compact
                      << " -> globalBack=" << globalBack << std::endl;
            assert(false);
        }
    }

    std::cout << "  OK!" << std::endl;
}

int DebugGlobalToCompactIdx(){
    // 1D: 10 bins
    TH1D h1("h1", "h1", 11, 0, 1);
    testHistogram(h1, "TH1D");

    // 2D: 5x7 bins
    TH2D h2("h2", "h2", 15, 0, 1, 7, 0, 1);
    testHistogram(h2, "TH2D");

    // 3D: 3x4x2 bins
    TH3D h3("h3", "h3", 3, 0, 1, 4, 0, 1, 7, 0, 1);
    testHistogram(h3, "TH3D");

    std::cout << "All tests passed!" << std::endl;
    return 0;
}
// All tests passed!