#ifndef ABSTRACTLLGSYSTEMSOLVERDEF
#define ABSTRACTLLGSYSTEMSOLVERDEF

#include <Eigen/Dense>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>
#include <math.h>
#include <vector>
using namespace std;
using namespace Eigen;

class AbstractLLGSystemSolver
{
private:
	double timeStepSize;
	double RStepSize;
    double ZStepSize;
	double initialTime;
	double finalTime;
	double magneticField;
	double anisotropyConstant;
	double dampingConstant;
	double DMconstant;
	double matR;
    double extR;
    double maxR;
    double matZ;
    double extZ;
    double maxZ;
    double minZ;
    double matMinZ;
    double matMaxZ;
    int Nr_max_mat;
    int Nr_max;
    int Nz_max;
    int Jz_min_mat;
    int Jz_in_mat;
    int Jz_max_mat;
    
public:
	void SetTimeStepSize(double dt);
	void SetRStepSize(double dr);
    void SetZStepSize(double dz);
	void SetTimeInterval(double t0, double t1);
    void SetRMax(double r_max);
    void SetZMin(double z_min);
    void SetZMax(double z_max);
    void SetZMatMin(double z_min_mat);
    void SetZMatMax(double z_max_mat);
	void SetRMatMax(double r_mat);
    void SetExteriorRLength();
    void SetMaterialZLength();
    void SetExteriorZLength();
    void SetVariousLengths();
	void SetMagneticField(double H);
	void SetAnisotropyConstant(double kappa);
	void SetDampingConstant(double alpha);
	void SetDMconstant(double lambda);
    void SetNr_max();
    void SetNz_max();
    void SetNr_max_mat();
    void SetJz_min_mat();
    void SetJz_in_mat();
    void SetJz_max_mat();
    void SetSpaceDiscretization();
    
	virtual void SolveEquation() = 0;
    int LI(int i, int j, int J);

	double GetInitialTime();
	double GetFinalTime();
	double GetTimeStepSize();
	double GetRStepSize();
    double GetZStepSize();
    double GetRLengthOfMaterial();
    double GetRExteriorLength();
    double GetRMax();
    double GetZLengthOfMaterial();
    double GetZExteriorLength();
    double GetZMax();
    double GetZMin();
    double GetZMatMax();
    double GetZMatMin();
    int GetNr_max();
    int GetNr_max_mat();
    int GetNz_max();
    int GetJz_min_mat();
    int GetJz_in_mat();
    int GetJz_max_mat();
	double GetMagneticField();
	double GetAnisotropyConstant();
	double GetDampingConstant();
	double GetDMconstant();
};

#endif
