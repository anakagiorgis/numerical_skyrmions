#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <math.h>

#include "AbstractLLGSystemSolver.h"
#include "LLG.h"

using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
	LLG Eq0;
    Eq0.SetTimeStepSize(1e-3);
    Eq0.SetRStepSize(0.1);
    Eq0.SetZStepSize(0.1);
    Eq0.SetRMax(100.0);
    Eq0.SetRMatMax(10.0);
    Eq0.SetZMin(-10.0);
    Eq0.SetZMax(12.0);
    Eq0.SetZMatMin(0.0);
    Eq0.SetZMatMax(2.0);
    Eq0.SetVariousLengths();
    Eq0.SetAnisotropyConstant(3.0);
    Eq0.SetMagneticField(0.0);
    Eq0.SetDampingConstant(0.5);
    Eq0.SetDMconstant(1.0);
    Eq0.SolveEquation();
    
	return 0;
}
