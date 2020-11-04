// Solve Optimal control problem for tumor diameter, T, number of infilitrating lympocytes, I,
// and concentration of of a chemotherapy drug, I.
// Use the Foward_backward sweep method - RK4 to solve the state and BRK4 for the adjoint
// Obj function J(u) = \int_{0}^{t_{end}} = (T+k*u_1^2)dt 
#include<iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <algorithm>
using namespace std;

int main()
{
    // Time and step size parameters 

    double tend = 24; // total time to run for (hrs) 
    int tsteps = 1000, itt = 1; // number of time steps, itteration count
    double h = tend / tsteps;
    double h2 = h / 2;

    // Initialise vectors to store solutions \\ 

    vector<double> T(tsteps + 1); // Vector of tumor diameter at all times 
    vector<double> I(tsteps + 1); // Vector of infiltrating lymphocytes
    vector<double> X(tsteps + 1); // Concentration of drug 
    vector<double> L1(tsteps + 1); // Vector of Lambda1 solutions
    vector<double> L2(tsteps + 1); // Vector of Lambda2 solutions
    vector<double> L3(tsteps + 1); // Vector of Lambda3 solutions
    vector<double> u(tsteps + 1); // Vector of u solutions

    vector<double> OT(tsteps + 1); // Store old solutions for 
    vector<double> OI(tsteps + 1); // Vector of infiltrating lymphocytes
    vector<double> OX(tsteps + 1); // Concentration of drug 
    vector<double> OL1(tsteps + 1); // Vector of Lambda1 solutions
    vector<double> OL2(tsteps + 1); // Vector of Lambda2 solutions
    vector<double> OL3(tsteps + 1); // Vector of Lambda3 solutions
    vector<double> Ou(tsteps + 1); // Vector of Lambda3 solutions
    vector<double> temp(tsteps + 1, 0.0); // Vector to copare i+1 output with i'th output for u
    vector<double> u1(tsteps + 1); // To be used to in calculation of u each step 

    T[0] = 25; // ICS
    I[0] = 200; // ICS
    X[0] = 25; // ICS

    // Initiate and define parameters and variables \\

    double R1 = 3.31E-2, R3 = 0.0001, s = 0.09, e = 0.09, mu1 = 0.125, mu2 = 0.7, z1 = 0.8,
        z2 = 0.1, beta = 0.0001, Ke = 0.63/12, delta = 0.1, kapa = 10e4, test = -1;
    double m11, m12, m13, m21, m22, m23, m31, m32, m33, m41, m42, m43;
    double n11, n12, n13, n21, n22, n23, n31, n32, n33, n41, n42, n43;
    double  temp1, temp2, temp3, temp4, temp5, temp6, temp7;
    double sumu = 0, sumT = 0, sumX = 0, sumI = 0, sumL1 = 0, sumL2 = 0, sumL3 = 0, sumOu = 0,
        sumOT = 0, sumOX = 0, sumOI = 0, sumOL1 = 0, sumOL2 = 0, sumOL3 = 0;

    // Main bit we run through 
    while (test < 0)
    {
        OT = T;
        OI = I;
        OX = X;
        OL1 = L1;
        OL2 = L2;
        OL3 = L3;
        Ou = u;

        for (int i = 0; i < tsteps; i++) //loop through time for state ODE 
        {
            // Forward Runge Kutta 4th order sceme for T,I,X
            m11 = R1 * T[i] - e * T[i] * I[i] / (I[i] + s * T[i]) - mu1 * T[i] * (1 + z1 * X[i] / (1 + X[i]));
            m12 = beta * T[i] + R3 * I[i] * (e * T[i] * I[i] / (I[i] + s * T[i])) - mu2 * I[i] * (1 + z2 * X[i] / (1 + X[i]));
            m13 = u[i] - Ke * X[i];


            m21 = R1 * (T[i] + h2 * m11) - e * ((T[i] + h2 * m11) * (I[i] + h2 * m12) / ((I[i] + h2 * m12) + (s * (T[i] + h2 * m11))) - mu1 * (T[i] + h2 * m11) * (1 + z1 * (X[i] + h2 * m13) / (1 + (X[i] + h2 * m13))));
            m22 = beta * (T[i] + h2 * m11) + R3 * (I[i] + h2 * m12) * (e * (T[i] + h2 * m11) * (I[i] + h2 * m12) / ((I[i] + h2 * m12) + s * (T[i] + h2 * m11))) - mu2 * (I[i] + h2 * m12) * (1 + z2 * (X[i] + h2 * m13) / (1 + (X[i] + h2 * m13)));
            m23 = 0.5 * (u[i] + u[i + 1]) - Ke * (X[i] + h2 * m13);


            m31 = R1 * (T[i] + h2 * m21) - e * ((T[i] + h2 * m21) * (I[i] + h2 * m22) / ((I[i] + h2 * m22) + (s * (T[i] + h2 * m21))) - mu1 * (T[i] + h2 * m21) * (1 + z1 * (X[i] + h2 * m23) / (1 + (X[i] + h2 * m23))));
            m32 = beta * (T[i] + h2 * m21) + R3 * (I[i] + h2 * m22) * (e * (T[i] + h2 * m21) * (I[i] + h2 * m22) / ((I[i] + h2 * m22) + s * (T[i] + h2 * m21))) - mu2 * (I[i] + h2 * m22) * (1 + z2 * (X[i] + h2 * m23) / (1 + (X[i] + h2 * m23)));
            m33 = 0.5 * (u[i] + u[i + 1]) - Ke * (X[i] + h2 * m23);

            m41 = R1 * (T[i] + h * m31) - e * ((T[i] + h * m31) * (I[i] + h * m32) / ((I[i] + h * m32) + (s * (T[i] + h * m31))) - mu1 * (T[i] + h * m31) * (1 + z1 * (X[i] + h * m33) / (1 + (X[i] + h * m33))));
            m42 = beta * (T[i] + h * m31) + R3 * (I[i] + h * m32) * (e * (T[i] + h * m31) * (I[i] + h * m32) / ((I[i] + h * m32) + s * (T[i] + h * m31))) - mu2 * (I[i] + h * m32) * (1 + z2 * (X[i] + h * m33) / (1 + (X[i] + h * m33)));
            m43 = u[i + 1] - Ke * (X[i] + h * m33);


            T[i + 1] = T[i] + (h / 6) * (m11 + 2 * m21 + 2 * m31 + m41);
            I[i + 1] = I[i] + (h / 6) * (m12 + 2 * m22 + 2 * m32 + m42);
            X[i + 1] = X[i] + (h / 6) * (m13 + 2 * m23 + 2 * m33 + m43);


        }


        // Calculate the error for u, T, I, X 
        for (int ii = 0; ii < tsteps; ii++) //loop through time for Adjoint ODE 
        {
            int  j = tsteps + 1 - ii; // for backward RK4 sceme 

            n11 = -(1 + L1[j]) * (R1 - e * (I[j] - s)) / pow((I[j] + s * T[j]), 2) - mu1 * (1 + (z1 * X[j]) / (1 + X[j])) + L2[j] * (beta + R3 * (I[j] - s) / pow((I[j] + s * T[j]), 2));
            n12 = -(L1[j] - e * ((T[j] - 1) / pow(I[j] + s * T[j], 2)) + L2[j] * R3 * e * T[j] * I[j] / (I[j] + s * T[j]) + R3 * I[j] * e * (T[j] - 1) / pow(I[j] + s * T[j], 2) - mu2 * (1 + z2 * X[j] / (1 + X[j])));
            n13 = -(L1[j] * (-(mu1 * T[j] * z1 - 1) / pow(1 + X[j], 2)) + L2[j] * (-(mu2 * I[j] * z2) / pow(1 + X[j], 2)) - L3[j] * (Ke));

            n21 = -(1 + (L1[j] - h2 * n11) * (R1 - e * ((0.5 * (I[j] + I[j - 1]) - s) / (0.5 * (I[j] + I[j - 1]) + s * (0.5 * pow(T[j] + T[j - 1], 2)))) - mu1 * ((1 + z1 * 0.5 * (X[j] + X[j - 1])) /
                (1 + 0.5 * (X[j] + X[j - 1])))) + (L2[j] - h2 * n12) * (beta + R3 * (0.5 * (I[j] + I[j - 1]) - s) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * pow(T[j] + T[j - 1], 2))));
            n22 = -((L1[j] - h2 * n11) - (e * (0.5 * (T[j] + T[j - 1]) - 1) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * pow(T[j] + T[j - 1], 2))) + (L2[j] - h2 * n12) * (R3 * e * 0.5 * (T[j] + T[j - 1])
                * 0.5 * (I[j] + I[j - 1]) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * (T[j] + T[j - 1])) + R3 * 0.5 * (I[j] + I[j - 1]) * e * (0.5 * (T[j] + T[j - 1]) - 1) / (0.5 * (I[j] + I[j - 1]) + s
                    * 0.5 * pow(T[j] + T[j - 1], 2)) - mu2 * (1 + z2 * 0.5 * (X[j] + X[j - 1]) / (1 + 0.5 * (X[j] + X[j - 1])))));
            n23 = -((L1[j] - h2 * n11) * (-1 * (mu1 * 0.5 * (T[j] + T[j - 1]) * z1 - 1) / (1 + 0.5 * pow(X[j] + X[j - 1], 2))) + (L2[j] - h2 * n12) * (-(mu2 * 0.5 * (I[j] + I[j - 1]) * z2 - 1) /
                (1 + 0.5 * pow(X[j] + X[j - 1], 2))) - (L3[j] - h2 * n13) * (Ke));


            n31 = -(1 + (L1[j] - h2 * n21) * (R1 - e * ((0.5 * (I[j] + I[j - 1]) - s) / (0.5 * (I[j] + I[j - 1]) + s * (0.5 * pow(T[j] + T[j - 1], 2)))) - mu1 * ((1 + z1 * 0.5 * (X[j] + X[j - 1]))
                / (1 + 0.5 * (X[j] + X[j - 1])))) + (L2[j] - h2 * n22) * (beta + R3 * (0.5 * (I[j] + I[j - 1]) - s) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * pow(T[j] + T[j - 1], 2))));
            n32 = -((L1[j] - h2 * n21) - (e * (0.5 * (T[j] + T[j - 1]) - 1) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * pow(T[j] + T[j - 1], 2))) + (L2[j] - h2 * n22) * (R3 * e * 0.5 * (T[j] + T[j - 1])
                * 0.5 * (I[j] + I[j - 1]) / (0.5 * (I[j] + I[j - 1]) + s * 0.5 * (T[j] + T[j - 1])) + R3 * 0.5 * (I[j] + I[j - 1]) * e * (0.5 * (T[j] + T[j - 1]) - 1) / (0.5 * (I[j] + I[j - 1])
                    + s * 0.5 * pow(T[j] + T[j - 1], 2)) - mu2 * (1 + z2 * 0.5 * (X[j] + X[j - 1]) / (1 + 0.5 * (X[j] + X[j - 1])))));
            n33 = -((L1[j] - h2 * n21) * (-(mu1 * 0.5 * (T[j] + T[j - 1]) * z1 - 1) / (1 + 0.5 * pow(X[j] + X[j - 1], 2))) + (L2[j] - h2 * n22) * (-(mu2 * 0.5 * (I[j] + I[j - 1]) * z2 - 1)
                / (1 + 0.5 * pow(X[j] + X[j - 1], 2))) - (L3[j] - h2 * n23) * (Ke));


            n41 = -(1 + (L1[j] - h * n31) * (R1 - e * ((I[j - 1] - s) / pow(I[j - 1] + s * T[j - 1], 2)) - mu1 * (1 + z1 * X[j - 1] / (1 + X[j - 1])))) + (L2[j] - h * n32) * (beta + R3 * (I[j - 1] - s)
                / pow(I[j - 1] + s * T[j - 1], 2));
            n42 = -((L1[j] - h * n31) - e * ((T[j - 1] - 1) / pow(I[j - 1] + s * T[j - 1], 2)) + (L2[j] - h * n32) * (R3 * e * T[j - 1] * I[j - 1] / (I[j - 1] + s * T[j - 1]) + R3 * I[j - 1] * e *
                (T[j - 1] - 1) / pow(I[j - 1] + s * T[j - 1], 2) - mu2 * (1 + z2 * X[j - 1] / (1 + X[j - 1]))));
            n43 = -((L1[j] - h * n31) * (-(mu1 * T[j - 1] * z1 - 1) / pow(1 + X[j - 1], 2)) + (L2[j] - h * n32) * (-(mu2 * I[j - 1] * z2 - 1) / pow(1 + X[j - 1], 2)) - (L3[j] - h * n33) * (Ke));

            L1[j - 1] = L1[j] - (h / 6) * (n11 + 2 * n21 + 2 * n31 + n41);
            L2[j - 1] = L2[j] - (h / 6) * (n12 + 2 * n22 + 2 * n32 + n42);
            L3[j - 1] = L3[j] - (h / 6) * (n13 + 2 * n23 + 2 * n33 + n43);
        }

        for (int j = 0; j < tsteps; j++) //Calculate 'u' for this itteration
        {
            temp[j] = -1 * L3[j] / (kapa);

            for (int jj = 0; jj < 7; jj++) // find min val in temp array 
            {
                if (u1[jj] > 50)
                    u1[jj] = 50;
                else
                    u1[jj] = 50; 
            }
            u[j] = 0.5 * (u[j] + Ou[j]);
        }

        // Calculate the error for u, T, I, X 
        for (int iii = 0; iii < tsteps; iii++)
        {
            sumu += abs(u[iii]);
            sumT += abs(T[iii]);
            sumX += abs(X[iii]);
            sumI += abs(I[iii]);
            sumL1 += abs(L1[iii]);
            sumL2 += abs(L2[iii]);
            sumL3 += abs(L3[iii]);
            sumOu += abs(Ou[iii]);
            sumOT += abs(OT[iii]);
            sumOX += abs(OX[iii]);
            sumOI += abs(OI[iii]);
            sumOL1 += abs(OL1[iii]);
            sumOL2 += abs(OL2[iii]);
            sumOL3 += abs(OL3[iii]);

        }
        

        temp1 = delta * sumu - sumOu;
        temp2 = delta * sumT - sumOT;
        temp3 = delta * sumX - sumOX;
        temp4 = delta * sumI - sumOI;
        temp5 = delta * sumL1 - sumOL1;
        temp6 = delta * sumL2 - sumOL2;
        temp7 = delta * sumL3 - sumOL3;

        double temparr[] = { temp1, temp2, temp3, temp4, temp5, temp6, temp7 };
        double min = 1000000;

        for (int j = 0; j < 7; j++) // find min val in temp array 
        {
            if (min > temparr[j])
                min = temparr[j];
            else
                min = min;
        }
        test = min;

        cout << " Min val =  " << test << "\n";
        if (itt > 10) // ensure we don't get stuck in infinite while loop
        {
            test = 2;
        }
        cout << " Itteration number  " << itt << "\n";
        itt = itt + 1; // keep track of the number of itterations
        cout  << temp1 << "\n" << temp2 << "\n" << temp3 << "\n" << temp4 << "\n" << temp5 << "\n" << temp6 << "\n" << temp7;

   

        
    }
}

