#include"General.h"


VectorXd CreateRandomVector(int size)
{
	VectorXd res(size);
	srand(time(NULL));
	for (int i = 0; i < size; i++)res(i) = rand(-1, 1);
	return res;
}
MatrixXd CreateRandomMatrix(int rows, int cols)
{
	srand(time(NULL));
	MatrixXd res(rows, cols);
	int size = rows * cols;
		for (int row = 0; row < rows; row++)
		{
			for (int col = 0; col < cols; col++)
			{
				res(row, col) = rand(-1, 1);
			}
		}
	return res;
}
double rand(double left, double right)
{
	return left + (right - left) * double(rand()) / RAND_MAX;
}
void DoSvd(MatrixXd& A, MatrixXd&U, MatrixXd& V, VectorXd& S)
{
	JacobiSVD<MatrixXd>svd(A, ComputeThinU | ComputeThinV);
	U = svd.matrixU();
	V = svd.matrixV();
	S = svd.singularValues();
}
VectorXd sigma_1(VectorXd& S)
{
	int size = S.size();
	double cap = 0.05 * S(0);
	VectorXd res(size);
	for (int i = 0; i < size; i++)
	{
		if (S(i) < cap)res(i) = 0;
		else res(i) = 1. / S(i);
	}
	return res;
}
double CalcDiscrepancy(MatrixXd& A, VectorXd& B, MatrixXd& x)
{
	MatrixXd temp = A * x;
	temp -= B;
	int size = temp.rows();
	double res = 0;
	for (int i = 0; i < size; i++)
	{
		res += temp(i, 0) * temp(i, 0);
	}
	return res/* / size */;
}
double CalcDiscrepancy(MatrixXd& x, VectorXd& CX)
{
	VectorXd temp;
	double res = 0;
	temp = x - CX;
	int size = temp.rows();
	for (int i = 0; i < size; i++)
	{
		res += temp(i, 0) * temp(i, 0);
	}
	return res /*/ size*/;
}


VectorXd AddNoise(VectorXd& B, float NoiseLevel)
{
	srand(clock());
	VectorXd res = B;
	int size = res.rows();
	for (int i = 0; i < size; i++)
	{
		res(i) = res(i) * rand(1 - NoiseLevel, 1 + NoiseLevel);
	}
	return res;
}


ostream& PrintMatrix(ostream& stream, MatrixXd& A, const char* MatrixTitle)
{
	stream << "\n\nMatrix "<<MatrixTitle<<":\n\n";
	
	int rows = A.rows();
	int cols = A.cols();
	if ((rows == 0) || (cols == 0))
	{
		stream << "Matrix "<<MatrixTitle<<" is empty\n\n";
		return stream;
	}


	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols; col++)
		{
			stream << setprecision(6) << setw(11) << setfill(' ') << A(row, col);
		}
		stream << endl;
	}
	return stream;
}
ostream& PrintVector(ostream& stream, VectorXd& A, const char* VectorTitle)
{
	stream << "\n\nVector " << VectorTitle << ":\n\n";
	int size = A.size();
	if (size == 0)
	{
		stream << "Vector " << VectorTitle << " is empty\n";
		return stream;
	}
	for (int i = 0; i < size; i++)stream << setprecision(6) << setw(11) << setfill(' ') << A(i) << endl;
	stream << endl;
	return stream;
}


void SolveLinearEquationsSystem()
{
	int rows(3), cols(5);
	cout << "\n\nWrite size rows x cols\n";

	cin >> rows >> cols;

	auto A = CreateRandomMatrix(rows, cols);
	//MatrixXd A(rows, cols);
	//A << 1, 2, 2, 4, 3, 5, 6, 7, 8, 9, 3, 2;
	auto CertainX = CreateRandomVector(cols);
	//VectorXd CertainX(cols);
	//CertainX << 1, 1, 1, 1, 1;
	VectorXd B = A * CertainX;
	//B = AddNoise(B, 0.10);
	PrintMatrix(cout, A, "A");
	PrintVector(cout, B, "B");
	PrintVector(cout, CertainX, "Certain X");

	MatrixXd U, V;
	VectorXd S;
	DoSvd(A, U, V, S);
	PrintMatrix(cout, U, "U");
	PrintMatrix(cout, V, "V");
	PrintVector(cout, S, "S");

	cout << "\n\n Casuallity: " << S(0) / S(min(rows, cols) - 1);
	S = sigma_1(S);
	PrintVector(cout, S, "Inversed S");

	MatrixXd x = V * S.asDiagonal();
	x *= U.transpose();
	x = x * B;
	PrintMatrix(cout, x, "x");

	cout << "\n\nDisperancy(dB): " << CalcDiscrepancy(A, B, x) << endl;
	cout << "\n\nDisperancy(dX): " << CalcDiscrepancy(x, CertainX) << endl;
}