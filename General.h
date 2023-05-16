#pragma once



#include<eigen3/Eigen/SVD>
#include<iostream>
#include<iomanip>
using namespace Eigen;
using namespace std;


VectorXd CreateRandomVector(int size);
MatrixXd CreateRandomMatrix(int rows, int cols);
double rand(double left, double right);
void DoSvd(MatrixXd& A, MatrixXd& U, MatrixXd& V, VectorXd& S);
VectorXd sigma_1(VectorXd& S);

ostream& PrintMatrix(ostream& stream, MatrixXd& A, const char* MatrixTitle = "");
ostream& PrintVector(ostream& stream, VectorXd& A, const char* VectorTitle = "");


void SolveLinearEquationsSystem();

double CalcDiscrepancy(MatrixXd& A, VectorXd& B, MatrixXd& x);