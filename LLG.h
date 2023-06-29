#ifndef LANDAULIFSHITZGILBERTDEF
#define LANDAULIFSHITZGILBERTDEF

#include "AbstractLLGSystemSolver.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <math.h>
#include <vector>
using namespace std;
using namespace Eigen;

class LLG: public AbstractLLGSystemSolver
{
public:
    void SetRandZ();
    void RightHandSide(double m1[], double m2[], double m3[], double rhs1[], double rhs2[], double rhs3[], double r0[]);
	void SolveEquation();
    double Trapz(double m3[], double r0[]);
    void InitialConditions(double m1[], double m2[], double m3[], double r0[]);
    void NeumannBoundaryConditions(double m1[], double m2[], double m3[]);
    void WriteM3(double m3[]);
    void WriteInitialM3(double m3[]);
    void WriteInitialM2(double m2[]);
};

#endif