#include "AbstractLLGSystemSolver.h"
#include <cmath>
#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
using namespace std;
using namespace Eigen;

// Function to find the vector index from a matrix index

int AbstractLLGSystemSolver::LI(int i, int j, int J)
{
    int ind = j + i*J;
    return ind;
}


// Set time step size method

void AbstractLLGSystemSolver::SetTimeStepSize(double dt)
{
	timeStepSize = dt;
}


// Set radial distance step size method

void AbstractLLGSystemSolver::SetRStepSize(double dr)
{
	RStepSize = dr;
}


// Set z-direction distance step size method

void AbstractLLGSystemSolver::SetZStepSize(double dz)
{
    ZStepSize = dz;
}


// Set time interval method

void AbstractLLGSystemSolver::SetTimeInterval(double t0, double t1)
{
	initialTime = t0;
	finalTime = t1;
}

// Set the maximum value of rho

void AbstractLLGSystemSolver::SetRMax(double r_max)
{
    maxR = r_max;
}

// Set the length of the material in terms of rho

void AbstractLLGSystemSolver::SetRMatMax(double r_mat)
{
    matR = r_mat;
}

// Set the length outside the material in terms of rho

void AbstractLLGSystemSolver::SetExteriorRLength()
{
    double r_mat = AbstractLLGSystemSolver::GetRLengthOfMaterial();
    double r_max = AbstractLLGSystemSolver::GetRMax();
    double r_ext = r_max - r_mat;
    extR = r_ext;
}

// Set the length of the material in terms of z 

void AbstractLLGSystemSolver::SetMaterialZLength()
{
    double z_mat_min = AbstractLLGSystemSolver::GetZMatMin();
    double z_mat_max = AbstractLLGSystemSolver::GetZMatMax();
    double z_mat = z_mat_max - z_mat_min;
    matZ = z_mat;
}


// Set the length outside the material in terms of z 

void AbstractLLGSystemSolver::SetExteriorZLength()
{
    double z_max = AbstractLLGSystemSolver::GetZMax();
    double z_mat_max = AbstractLLGSystemSolver::GetZMatMax();
    double z_ext = z_max - z_mat_max;
    extZ = z_ext;
}


void AbstractLLGSystemSolver::SetVariousLengths()
{
    AbstractLLGSystemSolver::SetExteriorRLength();
    AbstractLLGSystemSolver::SetMaterialZLength();
    AbstractLLGSystemSolver::SetExteriorZLength();
}


// Set the maximum value of z

void AbstractLLGSystemSolver::SetZMax(double z_max)
{
    maxZ = z_max;
}


// Set the minimum value of z

void AbstractLLGSystemSolver::SetZMin(double z_min)
{
    minZ = z_min;
}

// Set the value where the material starts in the z-direction

void AbstractLLGSystemSolver::SetZMatMin(double z_min_mat)
{
    matMinZ = z_min_mat;
}

// Set the value where the material ends in the z-direction

void AbstractLLGSystemSolver::SetZMatMax(double z_max_mat)
{
    matMaxZ = z_max_mat;
}

// Set the value of the external magnetic field (only in z-direction)

void AbstractLLGSystemSolver::SetMagneticField(double H)
{
	magneticField = H;
}

// Set the value of the anisotropy constant kappa (length scale of the problem)

void AbstractLLGSystemSolver::SetAnisotropyConstant(double kappa)
{
	anisotropyConstant = kappa;
}

// Set the value of the Gilbert damping constant alpha

void AbstractLLGSystemSolver::SetDampingConstant(double alpha)
{
	dampingConstant = alpha;
}

// Set the value of the D-M constant lambda

void AbstractLLGSystemSolver::SetDMconstant(double lambda)
{
	DMconstant = lambda;
}


void AbstractLLGSystemSolver::SetNr_max_mat()
{
    double R_ext = AbstractLLGSystemSolver::GetRLengthOfMaterial();
    double dr = AbstractLLGSystemSolver::GetRStepSize();
    Nr_max_mat = int(round(R_ext/dr));
}

void AbstractLLGSystemSolver::SetNr_max()
{
    double R_max = AbstractLLGSystemSolver::GetRMax();
    double dr = AbstractLLGSystemSolver::GetRStepSize();
    Nr_max = int(round(R_max/dr));
}


void AbstractLLGSystemSolver::SetNz_max()
{
    double Z_max = AbstractLLGSystemSolver::GetZMax();
    double Z_min = AbstractLLGSystemSolver::GetZMin();
    double Z_tot = Z_max - Z_min;
    double dz = AbstractLLGSystemSolver::GetZStepSize();
    Nz_max = int(round(Z_tot/dz));
}


void AbstractLLGSystemSolver::SetJz_in_mat()
{
    double Z_in = AbstractLLGSystemSolver::GetZLengthOfMaterial();
    double dz = AbstractLLGSystemSolver::GetZStepSize();
    Jz_in_mat = int(round(Z_in/dz));
}


void AbstractLLGSystemSolver::SetJz_min_mat()
{
    double Z_min_mat = AbstractLLGSystemSolver::GetZExteriorLength();
    double dz = AbstractLLGSystemSolver::GetZStepSize();
    Jz_min_mat = int(round(Z_min_mat/dz));
}


void AbstractLLGSystemSolver::SetJz_max_mat()
{
    Jz_max_mat = Jz_min_mat + Jz_in_mat;
}

void AbstractLLGSystemSolver::SetSpaceDiscretization()
{
    AbstractLLGSystemSolver::SetNr_max();
    AbstractLLGSystemSolver::SetNz_max();
    AbstractLLGSystemSolver::SetNr_max_mat();
    AbstractLLGSystemSolver::SetJz_in_mat();
    AbstractLLGSystemSolver::SetJz_min_mat();
    AbstractLLGSystemSolver::SetJz_max_mat();
}


// Getters

// Get time step size method

double AbstractLLGSystemSolver::GetTimeStepSize()
{
	return timeStepSize;
}

// Get radial direction step size method

double AbstractLLGSystemSolver::GetRStepSize()
{
	return RStepSize;
}

// Get z-direction step size method

double AbstractLLGSystemSolver::GetZStepSize()
{
    return ZStepSize;
}

// Get initial time method

double AbstractLLGSystemSolver::GetInitialTime()
{
	return initialTime;
}

// Get final time method

double AbstractLLGSystemSolver::GetFinalTime()
{
	return finalTime;
}

// Get the length of the material in terms of rho

double AbstractLLGSystemSolver::GetRLengthOfMaterial()
{
	return matR;
}

// Get the maximum value of rho

double AbstractLLGSystemSolver::GetRMax()
{
    return maxR;
}

// Get the length outside the material in terms of rho

double AbstractLLGSystemSolver::GetRExteriorLength()
{
	return extR;
}

// Get the length of the material in terms of z

double AbstractLLGSystemSolver::GetZLengthOfMaterial()
{
    return matZ;
}

// Get the length outside (below=above) the material in terms of z

double AbstractLLGSystemSolver::GetZExteriorLength()
{
    return extZ;
}

// Get the maximum value of z 

double AbstractLLGSystemSolver::GetZMax()
{
    return maxZ;
}

// Get the minimum value of z

double AbstractLLGSystemSolver::GetZMin()
{
    return minZ;
}

// Get the value of z where the material starts

double AbstractLLGSystemSolver::GetZMatMin()
{
    return matMinZ;
}

double AbstractLLGSystemSolver::GetZMatMax()
{
    return matMaxZ;
}

// Get external magnetic field method

double AbstractLLGSystemSolver::GetMagneticField()
{
	return magneticField;
}

// Get anisotropy constant method

double AbstractLLGSystemSolver::GetAnisotropyConstant()
{
	return anisotropyConstant;
}

// Get damping constant method

double AbstractLLGSystemSolver::GetDampingConstant()
{
	return dampingConstant;
}

// Get D-M constant method

double AbstractLLGSystemSolver::GetDMconstant()
{
	return DMconstant;
}


int AbstractLLGSystemSolver::GetNr_max()
{
    return Nr_max;
}

int AbstractLLGSystemSolver::GetNr_max_mat()
{
    return Nr_max_mat;
}

int AbstractLLGSystemSolver::GetJz_in_mat()
{
    return Jz_in_mat;
}

int AbstractLLGSystemSolver::GetJz_min_mat()
{
    return Jz_min_mat;
}

int AbstractLLGSystemSolver::GetJz_max_mat()
{
    return Jz_max_mat;
}
