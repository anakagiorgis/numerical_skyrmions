#define _USE_MATH_DEFINES

#include "LLG.h"
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <math.h>
#include <Eigen/IterativeLinearSolvers>
#include <vector>
using namespace std;
using namespace Eigen;

// We solve the Landau-Lifshitz-Gilbert equation in three dimensions (x,y,z) with axial
// symmetry so the magnetisation is a function of axial distance rho and z so the problem
// reduces to two dimensions (rho,z) and m_1, m_2, m_3 correspond to m_rho, m_phi, m_z.
// We solve this system by an explicit scheme and then project the solution to the
// unit sphere to enforce the constraint.

// We also evaluate a double integral using the trapezium rule
// to evaluate an energy to determine when the system has converged to a static
// solution.


double dr = 0.1;
double dz = 0.1;
int const Nr = 100;
int const Nz = 20;
double m1[(Nr+1)*(Nz+1)];
double m2[(Nr+1)*(Nz+1)];
double m3[(Nr+1)*(Nz+1)];

double r0[Nr+1];
double z0[Nz+1];

void LLG::SetRandZ()
{
    double z_min = LLG::GetZMatMin();
    
    // Set the rho vector
    for (int jr=0; jr<=Nr; jr++)
    {
        r0[jr] = jr*dr;
    }
    // Set the z vector
    for (int jz=0; jz<=Nz; jz++)
    {
        z0[jz] = z_min + jz*dz;
    }
}

void LLG::InitialConditions(double m1[], double m2[], double m3[], double r0[])
{
    for (int jz=0; jz<=Nz; jz++)
    {
        m1[LLG::LI(0,jz,Nz+1)] = 0.0;
        m2[LLG::LI(0,jz,Nz+1)] = 0.0;
        m3[LLG::LI(0,jz,Nz+1)] = -1.0;
        for (int jr=1; jr<=Nr; jr++)
        {
            m1[LLG::LI(jr,jz,Nz+1)] = 0.0;
            m3[LLG::LI(jr,jz,Nz+1)] = 1.0 - 2.0*exp(-r0[jr]*r0[jr]);
            m2[LLG::LI(jr,jz,Nz+1)] = sqrt(1.0 - m3[LLG::LI(jr,jz,Nz+1)]*m3[LLG::LI(jr,jz,Nz+1)]);
        }
    }
}

void LLG::NeumannBoundaryConditions(double m1[], double m2[], double m3[])
{
    for (int jz=0; jz<=Nz; jz++)
    {
        m1[LLG::LI(0,jz,Nz+1)] = 0.0;
        m2[LLG::LI(0,jz,Nz+1)] = 0.0;
        m3[LLG::LI(0,jz,Nz+1)] = -1.0;
        
        m1[LLG::LI(Nr,jz,Nz+1)] = m1[LLG::LI(Nr-1,jz,Nz+1)];
        m2[LLG::LI(Nr,jz,Nz+1)] = m2[LLG::LI(Nr-1,jz,Nz+1)];
        m3[LLG::LI(Nr,jz,Nz+1)] = m3[LLG::LI(Nr-1,jz,Nz+1)];
    }
    
    for (int jr=1; jr<=Nr-1; jr++)
    {
        m1[LLG::LI(jr,0,Nz+1)] = m1[LLG::LI(jr,1,Nz+1)];
        m2[LLG::LI(jr,0,Nz+1)] = m2[LLG::LI(jr,1,Nz+1)];
        m3[LLG::LI(jr,0,Nz+1)] = m3[LLG::LI(jr,1,Nz+1)];
        
        m1[LLG::LI(jr,Nz,Nz+1)] = m1[LLG::LI(jr,Nz-1,Nz+1)];
        m2[LLG::LI(jr,Nz,Nz+1)] = m2[LLG::LI(jr,Nz-1,Nz+1)];
        m3[LLG::LI(jr,Nz,Nz+1)] = m3[LLG::LI(jr,Nz-1,Nz+1)];
    }
}


void LLG::WriteInitialM3(double m3[])
{
    ofstream initialfile;
    initialfile.open("/home/giorgisk/Documents/workspace/LLGarrays/Initial.txt");
    
    for (int jr=0; jr<=Nr; jr++)
    {
        for (int jz=0; jz<=Nz; jz++)
        {
            initialfile << m3[LLG::LI(jr,jz,Nz+1)] << " ";
        }
            
        initialfile << "\n";
    }
        
    initialfile.close();
}

void LLG::WriteInitialM2(double m2[])
{
    ofstream initialm2file;
    initialm2file.open("/home/giorgisk/Documents/workspace/LLGarrays/Initialm2.txt");
    
    for (int jr=0; jr<=Nr; jr++)
    {
        for (int jz=0; jz<=Nz; jz++)
        {
            initialm2file << m2[LLG::LI(jr,jz,Nz+1)] << " ";
        }
            
        initialm2file << "\n";
    }
        
    initialm2file.close();
}

void LLG::WriteM3(double m3[])
{
    
    // Create a file to store the values of m_3 (m_z)
    ofstream newfile;
    newfile.open("/home/giorgisk/Documents/workspace/LLGarrays/Landau.txt");
    
    // Store the final values of m_3 (m_z), after the loop has been exited and the algorithm has stopped. 
    for (int jr=0; jr<=Nr; jr++)
    {
        for (int jz=0; jz<=Nz; jz++)
        {
            newfile << m3[LLG::LI(jr,jz,Nz+1)] << " ";
        }
            
        newfile << "\n";
    }
        
    newfile.close();
}
 
 
double LLG::Trapz(double m3[], double r0[])
{
    
    double sum = 0.0;
    
    for (int jr=0; jr<=Nr; jr++)
    {
        for (int jz=0; jz<=Nz; jz++)
        {
            sum += (1 - m3[LLG::LI(jr,jz,Nz+1)])*r0[jr];
        }
    }
    double totalmag;
    totalmag = 2*M_PI*dr*dz*sum;
    return totalmag;
}

void LLG::RightHandSide(double m1[], double m2[], double m3[], double rhs1[], double rhs2[], double rhs3[], double r0[])
{
    
    // Create variables to store the Laplacian of m_1, m_2, m_3.
	double Lm_1;
	double Lm_2;
	double Lm_3;
    
    // Create variables to store the Curl of m in the rho, phi,
    // and z directions.
    double Curl_1;
    double Curl_2;
    double Curl_3;
    
    // Create variables to store the Effective Field f.
    double f1;
    double f2;
    double f3;
    
    // Create a variable to store the inner product of f and m. 
    double dotproduct;
    
    // Create variables to store the cross product of m and f.
    double crossproduct1;
    double crossproduct2;
    double crossproduct3;
    
    double d2rho = dr*dr;
    double d2z = dz*dz;
    
    // Get the values of the constants to be used.
    double kappa = LLG::GetAnisotropyConstant();
    
    double alpha = LLG::GetDampingConstant();
    double lambda = LLG::GetDMconstant();
    
    double H = LLG::GetMagneticField();
    
    for (int jz=1; jz<=Nz-1; jz++)
            {
                for (int jr=1; jr<=Nr-1; jr++)
                {
                    // Compute the values of the Laplacians.
                    Lm_1 
 = (m1[LLG::LI(jr+1,jz,Nz+1)] - 2.0*m1[LLG::LI(jr,jz,Nz+1)] + m1[LLG::LI(jr-1,jz,Nz+1)])/(d2rho)
 + (m1[LLG::LI(jr+1,jz,Nz+1)]- m1[LLG::LI(jr,jz,Nz+1)])/(dr*r0[jr])
 + (m1[LLG::LI(jr,jz+1,Nz+1)] - 2.0*m1[LLG::LI(jr,jz,Nz+1)] + m1[LLG::LI(jr,jz-1,Nz+1)])/(d2z)
 - (m1[LLG::LI(jr,jz,Nz+1)]/(r0[jr]*r0[jr]));
 
                    Lm_2
 = (m2[LLG::LI(jr+1,jz,Nz+1)] - 2.0*m2[LLG::LI(jr,jz,Nz+1)] + m2[LLG::LI(jr-1,jz,Nz+1)])/(d2rho)
 + (m2[LLG::LI(jr+1,jz,Nz+1)] - m2[LLG::LI(jr,jz,Nz+1)])/(dr*r0[jr])
 + (m2[LLG::LI(jr,jz+1,Nz+1)] - 2.0*m2[LLG::LI(jr,jz,Nz+1)] + m2[LLG::LI(jr,jz-1,Nz+1)])/(d2z)
 - (m2[LLG::LI(jr,jz,Nz+1)]/(r0[jr]*r0[jr]));
 
                    Lm_3
 = (m3[LLG::LI(jr+1,jz,Nz+1)] - 2.0*m3[LLG::LI(jr,jz,Nz+1)] + m3[LLG::LI(jr-1,jz,Nz+1)])/(d2rho)
 + (m3[LLG::LI(jr+1,jz,Nz+1)] - m3[LLG::LI(jr,jz,Nz+1)])/(dr*r0[jr])
 + (m3[LLG::LI(jr,jz+1,Nz+1)] - 2.0*m3[LLG::LI(jr,jz,Nz+1)] + m3[LLG::LI(jr,jz-1,Nz+1)])/(d2z);
                    
                    // Compute the values of the Curls.
                    // Curl_1 = - (m2[LLG::LI(jr,jz+1,Nz+1)] - m2[LLG::LI(jr,jz-1,Nz+1)])/(2.0*dz);
                    // Curl_2 = (m1[LLG::LI(jr,jz+1,Nz+1)] - m1[LLG::LI(jr,jz-1,Nz+1)])/(2.0*dz) - (m3[LLG::LI(jr+1,jz,Nz+1)] - m3[LLG::LI(jr-1,jz,Nz+1)])/(2.0*dr);
                    
                    Curl_2
 = - (m3[LLG::LI(jr+1,jz,Nz+1)] - m3[LLG::LI(jr,jz,Nz+1)])/(dr);
 
                    Curl_3
 = (m2[LLG::LI(jr+1,jz,Nz+1)] - m2[LLG::LI(jr,jz,Nz+1)])/(dr)
 + (m2[LLG::LI(jr,jz,Nz+1)]/(r0[jr]));
                    
                    // Compute the values of the Effective Field.
                    f1 = Lm_1;
                    f2 = Lm_2 - 2.0*lambda*Curl_2;
                    f3 = Lm_3 - 2.0*lambda*Curl_3 + kappa*m3[LLG::LI(jr,jz,Nz+1)] + H;
                    
                    // Compute the dot product of m with f, <m,f>.
                    dotproduct
 = m1[LLG::LI(jr,jz,Nz+1)]*f1 + m2[LLG::LI(jr,jz,Nz+1)]*f2 + m3[LLG::LI(jr,jz,Nz+1)]*f3;
                    
                    // Compute the cross product of -m and f, -mxf.
                    crossproduct1 = m3[LLG::LI(jr,jz,Nz+1)]*f2 - m2[LLG::LI(jr,jz,Nz+1)]*f3;
                    crossproduct2 = m1[LLG::LI(jr,jz,Nz+1)]*f3 - m3[LLG::LI(jr,jz,Nz+1)]*f1;
                    crossproduct3 = m2[LLG::LI(jr,jz,Nz+1)]*f1 - m1[LLG::LI(jr,jz,Nz+1)]*f2;
                    
                    // Compute the values of the RHS using the above calculations.
                    rhs1[LLG::LI(jr,jz,Nz+1)] = crossproduct1 + alpha*(f1 - dotproduct*m1[LLG::LI(jr,jz,Nz+1)]);
                    rhs2[LLG::LI(jr,jz,Nz+1)] = crossproduct2 + alpha*(f2 - dotproduct*m2[LLG::LI(jr,jz,Nz+1)]);
                    rhs3[LLG::LI(jr,jz,Nz+1)] = crossproduct3 + alpha*(f3 - dotproduct*m3[LLG::LI(jr,jz,Nz+1)]);
                    
                }
            }
    
}

void LLG::SolveEquation()
{
    
    
    // Get the values of the time step dt.
    double dt = LLG::GetTimeStepSize();
    
    // Set rho
    SetRandZ();
    
    
    // Test print
    cout << dt/(dr*dr*dz*dz) << "\n";
    
    LLG::InitialConditions(m1, m2, m3, r0);
    
    
    LLG::WriteInitialM3(m3);
    LLG::WriteInitialM2(m2);
    
    // test print
    // cout << "This point has been reached" << "\n";
    
    // Create a variable to store the norm of the magnetisation m, the square root
	// of the sum of the squares of m_1, m_2, m_3.
    double norm;
    
    // Create a vector to store the values of the two integrals we are going to compute.
	// The first integral is the integral of (1-m_3) before the algorithm has computed
	// the value of m_3, and the second integral is the integral of (1-m_3) after the
	// computation. 
	double I[2];
    
    // Initialise the integrals so their difference is bigger than the tolerance
	// (1e-8) and the criterion for stopping the algorithm is not satisfied.
	I[1] = 1;
	I[0] = 0;
    
    double rhs1[(Nr+1)*(Nz+1)];
    double rhs2[(Nr+1)*(Nz+1)];
    double rhs3[(Nr+1)*(Nz+1)];
    
    // test print
    int i=0;
    
    double time=0.0;
    
    double tinfo=2.0;
    
    double c;
    
    if (fmod(time,tinfo)<dt){
                cout << time << "\n";
        
            }

    // Main Loop
    
    while (time<1000.0){
        // Criterion to stop the algorithm

		// is that the energy from one time step to the next has changed by a small
		// amount therefore it is reasonable to expect that a static solution has been
		// reached.

		// For every time the criterion is not met, the values of m_1, m_2, m_3 are
		// updated for a next time step according to the explicit approximation of
		// the Landau-Lifshitz-Gilbert equation.
        
            
            
            // Impose boundary conditions 
            
            LLG::NeumannBoundaryConditions(m1, m2, m3);
            
            
            LLG::RightHandSide(m1, m2, m3, rhs1, rhs2, rhs3, r0);
            
            
            
            // Loop to evaluate the Landau-Lifshitz-Gilbert equation for one time step
            // for all the values of rho and z.
            
            // Euler Method 
            
            for (int jz=1; jz<=Nz-1; jz++)
            {
                for (int jr=1; jr<=Nr-1; jr++)
                {
                    
                
                    
                    // Use the above values to calculate m_rho, m_phi and m_z. 
                    m1[LLG::LI(jr,jz,Nz+1)]
 = m1[LLG::LI(jr,jz,Nz+1)] + dt*rhs1[LLG::LI(jr,jz,Nz+1)];
 
                    m2[LLG::LI(jr,jz,Nz+1)]
 = m2[LLG::LI(jr,jz,Nz+1)] + dt*rhs2[LLG::LI(jr,jz,Nz+1)];
 
                    m3[LLG::LI(jr,jz,Nz+1)]
 = m3[LLG::LI(jr,jz,Nz+1)] + dt*rhs3[LLG::LI(jr,jz,Nz+1)];
                    
                     
                    
                    // Project the solution to the unit sphere to enforce the constraint.
                    norm
 = sqrt(m1[LLG::LI(jr,jz,Nz+1)]*m1[LLG::LI(jr,jz,Nz+1)]
 + m2[LLG::LI(jr,jz,Nz+1)]*m2[LLG::LI(jr,jz,Nz+1)]
 + m3[LLG::LI(jr,jz,Nz+1)]*m3[LLG::LI(jr,jz,Nz+1)]);
                    
                    m1[LLG::LI(jr,jz,Nz+1)] = m1[LLG::LI(jr,jz,Nz+1)]/norm;
                    m2[LLG::LI(jr,jz,Nz+1)] = m2[LLG::LI(jr,jz,Nz+1)]/norm;
                    m3[LLG::LI(jr,jz,Nz+1)] = m3[LLG::LI(jr,jz,Nz+1)]/norm;
                    
                    
                }
            }
            
            time = time + dt;
            
            
            
            // Criterion + prints
            if (fmod(time,tinfo)<dt){
                I[1] = LLG::Trapz(m3, r0);
                c = abs(I[1] - I[0])/tinfo;
                cout << time << " " << I[1] << " " << c << "\n";
                if (c<1e-8){
                    break;
                }
                I[0] = I[1];
        
            }
            
            
        i+=1;
        }
        
        // Store the final values of m_3 (m_z).
        LLG::WriteM3(m3);
    
}