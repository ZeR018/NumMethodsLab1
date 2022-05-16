#pragma once
#include "MMN.h"
#include <iostream>

MMN::MMN()
{
	pi = 3.14;
}


MMN::MMN(int N, int M, TVector<double> XBorder, TVector<double> YBorder)
{
	pi = 2 * asin(1);
	n = N; m = M;

	xBorder = XBorder; yBorder = YBorder;

	h = (xBorder[xBorder.Size() - 1] - xBorder[0]) / (1.0 * n);
	k = (yBorder[yBorder.Size() - 1] - yBorder[0]) / (1.0 * m);
	hE = -1 / pow(h, 2);
	kE = -1 / pow(k, 2);
	A = -2 * (hE + kE);

	V = TMatrix<double>(m + 1);
	U = TMatrix<double>(m + 1);
	for (int j = 0; j <= m; j++)
	{
		V[j] = TVector<double>(n + 1);
		U[j] = TVector<double>(n + 1);
		for (int i = 0; i <= n; i++)
		{
			V[j][i] = 0.0;
			U[j][i] = 0.0;
		}
	}
	R = V;
	F = V;
	FunctionInicialisation();
	Inicialisation();
}


double MMN::F_Function(double x, double y)
{
	double res;

	res = - 4*exp(1 - pow(x,2) - pow(y, 2)) * (1 - pow(x, 2) - pow(y, 2));

	return res;
}


void MMN::FunctionInicialisation()
{
	double x, y = yBorder[0] + k;
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			F[j][i] = F_Function(x, y);
			x += h;
		}
		y += k;
	}
}


double MMN::ExactSolution(double x, double y)
{
	double res;

	res = exp(1 - pow(x, 2) - pow(y, 2));

	return res;
}


double MMN::XInicialConditions(double x, int Num)
{
	double res = 0.0;
	switch (Num)
	{
	case 1:
		res = exp(-pow(x,2));
		break;
	case 2:
		res = exp(-pow(x, 2));
		break;
	}
	return res;
}


double MMN::YInicialConditions(double y, int Num)
{
	double res = 0.0;
	switch (Num)
	{
	case 1:
		res = exp(-pow(y, 2));
		break;
	case 2:
		res = exp(-pow(y, 2));
		break;
	}
	return res;
}


void MMN::Inicialisation()
{
	double x, y = yBorder[0];
	for (int j = 0; j <= m; j++)
	{
		V[j][0] = YInicialConditions(y, 1);
		V[j][n] = YInicialConditions(y, 2);
		U[j][0] = YInicialConditions(y, 1);
		U[j][n] = YInicialConditions(y, 2);
		y += k;
	}
	x = xBorder[0] + h;
	for (int i = 1; i < n; i++)
	{
		V[0][i] = XInicialConditions(x, 1);
		V[m][i] = XInicialConditions(x, 2);
		U[0][i] = XInicialConditions(x, 1);
		U[m][i] = XInicialConditions(x, 2);
		x += h;
	}
}


TVector<double> MMN::MethodAccuracy(TVector<double> eps, TVector<int> MaxIterations, char Name) //�������� �������
{
	return 0;
}


TVector<double> MMN::MethodError(double eps, int MaxIterations) //����������� �������
{
	return 0;
}


void MMN::VectorNevyazki()
{
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i - 1] + V[j][i + 1]) + kE * (V[j - 1][i] + V[j + 1][i]) + F[j][i];
		}
}

double MMN::NevyazkaInf()
{
	double res = 0, temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = fabs(R[j][i]);
			if (temp > res)
				res = temp;
		}

	return res;
}

double MMN::NevyazkaEvkl()
{
	double res = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			res += pow(R[j][i], 2);
		}

	return sqrt(res);
}

void MMN::XInterpolation()
{
	double x, a, b;
	a = xBorder[0]; b = xBorder[xBorder.Size() - 1];
	for (int j = 1; j < m; j++)
	{
		x = xBorder[0] + h;
		for (int i = 1; i < n; i++)
		{
			V[j][i] = ((x - a) * V[j][n] - (x - b) * V[j][0]) / (b - a);
			x += h;
		}
	}
}

void MMN::YInterpolation()
{
	double y, c, d;
	c = yBorder[0]; d = yBorder[yBorder.Size() - 1];
	for (int i = 1; i < n; i++)
	{
		y = yBorder[0] + k;
		for (int j = 1; j < m; j++)
		{
			V[j][i] = ((y - c) * V[m][i] - (y - d) * V[0][i]) / (d - c);
			y += k;
		}
	}
}


void MMN::SaveGrid(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int i = 0; i < xBorder.Size(); i++)
		file << xBorder[i] << endl;
	file << h << endl;

	for (int i = 0; i < yBorder.Size(); i++)
		file << yBorder[i] << endl;
	file << k << endl;

	file.close();

}

void MMN::SaveData(string s)
{
	ofstream file(s);

	file << n << endl << m << endl;

	for (int j = 0; j <= m; j++)
	{

		for (int i = 0; i <= n; i++)
		{
			file << V[j][i] << endl;//"\t";
		}
		//std::cout << endl;
	}
	file.close();

}

void MMN::SaveExData(string s)
{
	ofstream file(s);
	double res;

	for (int j = 0; j <= m; j++)
	{
		for (int i = 0; i <= n; i++)
		{
			file << U[j][i] << endl;//"\t";
		}
		//std::cout << endl;
	}
	file.close();
}
//------------------------------------------------------------

void MMN::SetParams()
{
	double Arr = 0, ArAr = 0;
	double temp;
	for (int j = 1; j < m; j++)
	{
		for (int i = 1; i < n; i++)
		{
			temp = A * R[j][i] + hE * (R[j][i + 1] + R[j][i - 1]) + kE * (R[j + 1][i] + R[j - 1][i]);
			Arr += temp * R[j][i];
			ArAr += pow(temp, 2);
		}
	}
	tau = Arr / ArAr;
}

double MMN::Runner()
{
	double temp;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			R[j][i] = A * V[j][i] + hE * (V[j][i + 1] + V[j][i - 1]) + kE * (V[j + 1][i] + V[j - 1][i]) + F[j][i];
		}
			
	SetParams();
	double accurancy = 0;
	for (int j = 1; j < m; j++)
		for (int i = 1; i < n; i++)
		{
			temp = V[j][i];
			//
			V[j][i] = V[j][i] - tau * R[j][i];
			
			//
			temp = fabs(V[j][i] - temp);
			if (temp > accurancy)
				accurancy = temp;
		}
	
	return accurancy;
}

TVector<double> MMN::SolvingTestTask(TVector<double> eps, TVector<int> MaxIterations, char Name)
{
	// ��� ����������� ����
	/*MMN Solution(2 * n, 2 * m, xBorder, yBorder);
	switch (Name)
	{
	case 'X':
		Solution.XInterpolation();
		break;
	case 'Y':
		Solution.YInterpolation();
		break;
	}*/

	TVector<double> accurancy(2); //�������� ������
	//����� ����� � ����
	for (int i = 0; i < 2; i++)
		accurancy[i] = 1 + eps[i];

	TVector<int> IterationsCount(2);
	//
	while ((accurancy[0] > eps[0]) && (IterationsCount[0] < MaxIterations[0]))
	{
		accurancy[0] = Runner();
		IterationsCount[0]++;
	}

	// ��� ����������� ����
	/*while ((accurancy[1] > eps[1]) && (IterationsCount[1] < MaxIterations[1]))
	{
		accurancy[1] = Solution.Runner();
		IterationsCount[1]++;
	}*/
	//
	

	// ������ �������
	double y0 = yBorder[0];
	for (int j = 1; j < m; j++)
	{
		double x0 = xBorder[0];
		y0 += k;
		for (int i = 1; i < n; i++)
		{
			x0 += h;
			U[j][i] = ExactSolution(x0, y0);
		}
	}

	SaveData("MainSolutA.txt"); //��������� ���� ��������� ��������
	//Solution.SaveData("SupSolutA.txt");  //��������� ���� � ��������� �������� �� ��������������� �����
	SaveGrid("DifferenceA.txt"); //��������� ���� ��������� �����������
	SaveExData("ExData.txt"); //��������� ���� � �������� ��������

	ofstream Difference("DifferenceA.txt", ios::app); //������ � ���� ����� ���������� � ���������� ��������

	double error = 0; //�������� �������
	double sup; //��������������� ����������
	int ix, jy;
	for (int j = 0; j <= m; j++)
		for (int i = 0; i <= n; i++)
		{
			sup = fabs(V[j][i] - U[j][i]);
			Difference << sup << endl;
			if (sup > error)
			{
				error = sup;
				ix = i;
				jy = j;
			}
		}

	double temp;
	VectorNevyazki();
	//Solution.VectorNevyazki();

	sup = NevyazkaEvkl();
	//temp = Solution.NevyazkaEvkl();
	temp = 0;
	TVector<double> result(9);
	result[0] = accurancy[0];			//�������� ������ �� ����� (n+1,m+1)
	result[1] = IterationsCount[0];		//���������� �������� �� ����� (n+1,m+1)
	result[2] = accurancy[1];			//�������� ������ �� ����� (2n+1,2m+1)
	result[3] = IterationsCount[1];		//���������� �������� �� ����� (2n+1,2m+1)
	result[4] = error;					//�������� �������
	result[5] = sup;					//��������� ����� ������� �� �������� �����
	result[6] = temp;					//��������� ����� ������� �� ��������������� �����
	result[7] = xBorder[0] + ix * h;	//�������� x � ����� ������ �����
	result[8] = yBorder[0] + jy * k;	//�������� y � ����� ������ �����

	/*std::cout << "���������� ��������: " << result[1] << "\n";
	std::cout << "�������� ������: " << result[0] << "\n";
	std::cout << "�������� �������: " << result[4] << "\n";
	std::cout << "��������� ����� �������: " << result[5] << "\n";*/

	cout << "################################################# ��������� ��������� #################################################" << endl;
	cout << endl;
	cout << "��� ������� ������ �� ��������� [" << xBorder[0] << ", " << xBorder[1] << "]" << endl;
	cout << "������������ ����������� ����� � ������ ��������� �� x n = " << n << " � ������ ��������� �� y m = " << m << "." << endl;
	cout << "�����: ����� ����������� �������, �������� ���: " << tau << " (������� �� ��������� ��������)." << endl;
	cout << "�������� ��������� �� �������� = " << eps[0] << " � �� ����� �������� = " << MaxIterations[0] << "." << endl;
	cout << endl;
	cout << "####################################################### ������� #######################################################" << endl << endl;
	cout << "���������� �������� S: " << result[1] << ",  ����������� ��������: " << result[0] << "." << endl;
	cout << "����� ������ � ��������: " << result[5] << " (��� ������� ���� ������������ ��������� �����)." << endl;
	cout << "����������� ������� �������� ������: " << result[4] << "." << endl;
	cout << "������������ ���������� ������� � ���������� ������� ����������� � ���� x = " << result[7] << " y = " << result[8] << endl;
	cout << "� ������� ���������� ����������� ������������ ";
	if (Name == 'X')
	{
		std::cout << "������������ �� X." << endl;
	}
	if (Name == 'Y')
	{
		std::cout << "������������ �� Y." << endl;
	}
	//std::cout << "�������� ��������� ��� ������� �� ������ ����������� �������� ������� ����." << endl;
	return result;
}

