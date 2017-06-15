#include<iostream>
#include<string>
#include<cmath>
#include<vector>
#include<cmath>
#include<locale>
using namespace std;
__forceinline double original_func(double x)
{
	double a = (13.0 + 3 * x*x) / (2 * x*x*x);
	return pow(a, 1.0 / 3.0);
}
__forceinline double func(double x, double y)
{
	return (1-x*y*y*y)/(x*x*y*y);
}
double delt(double x, double y,double h)
{
	double k1 = h*func(x, y);
	double k2 = h*func(x + h / 2.0, y + k1 / 2.0);
	double k3 = h*func(x + h / 2.0, y + k2 / 2.0);
	double k4 = h*func(x + h, y + k3);
	double del = 1.0 / 6 * (k1 + 2.0 * k2+2.0*k3 + k4);
	return del;
}
double func_proiz_2(double x, double y,double y_1)
{
	double rez = (x*y*y*y*y - 2 - func(x, y)*(2 * x + x*x*y*y*y))/(x*x*x*y*y*y);
	return rez;
}
double func_proiz_3(double x, double y, double y_1, double y_2)
{
	return (6.0 - 2 * x*y*y*y) / (x*x*x*x*y*y) + y_1*((8.0 + 2.0 * x*y*y*y) / (x*x*x*y*y*y)) + y_1*y_1 *( 6.0 / (x*x*y*y*y*y)) + y_2*((-2 - x*y*y*y) / (x*x*y*y*y));
}
//дл€ 2 степени
double formula_1(double h, double y, double y_1_i, double y_1_ip1, double y_2_i, double y_2_ip1)
{
	return y + h / 2 * (y_1_i + y_1_ip1) + h*h / 12 * (y_2_i - y_2_ip1);
}
double formula_2(double h, double y_i, double y_im1, double y_im2, double y_2_i,double y_2_im1)
{
	return y_im2 + 3.0 * (y_i - y_im1) + h*h*(y_2_i - y_2_im1);
}
decltype(auto) method_1(double x0, double y0, double h, double eps, int n)
{
	vector<vector<double>> otrezok(n);
	for (auto& x : otrezok)
	{
		x.resize(4);
	}
	otrezok[0][0] = x0;
	otrezok[0][1] = y0;
	otrezok[0][2] = func(otrezok[0][0], otrezok[0][1]);
	otrezok[0][3] = func_proiz_2(otrezok[0][0], otrezok[0][1], otrezok[0][2]);
	for (int i = 1; i < 3; i++)
	{
		otrezok[i][0] = otrezok[i - 1][0] + h;
		otrezok[i][1] = otrezok[i - 1][1] + delt(otrezok[i - 1][0], otrezok[i - 1][1], h);
		otrezok[i][2] = func(otrezok[i][0], otrezok[i][1]);
		otrezok[i][3] = func_proiz_2(otrezok[i][0], otrezok[i][1], otrezok[i][2]);
	}
	vector<double>priblijenie_1(4);
	vector<double>priblijenie_2(4);
	vector<double>priblijenie_3(4);
	for (int i = 2; i < n - 1; i++)
	{
		while (true)
		{
			//
			priblijenie_1[0] = otrezok[i][0] + h;
			priblijenie_1[1] = formula_2(h, otrezok[i][1], otrezok[i - 1][1], otrezok[i - 2][1], otrezok[i - 1][3], otrezok[i - 2][3]);
			priblijenie_1[2] = func(priblijenie_1[0], priblijenie_1[1]);
			priblijenie_1[3] = func_proiz_2(priblijenie_1[0], priblijenie_1[1], priblijenie_1[2]);

			//
			priblijenie_2[0] = otrezok[i][0] + h;
			priblijenie_2[1] = formula_1(h, otrezok[i][1], otrezok[i][2], priblijenie_1[2], otrezok[i][3], priblijenie_1[3]);
			priblijenie_2[2] = func(priblijenie_2[0], priblijenie_2[1]);
			priblijenie_2[3] = func_proiz_2(priblijenie_2[0], priblijenie_2[1], priblijenie_2[2]);

			//
			priblijenie_3[0] = otrezok[i][0] + h;
			priblijenie_3[1] = formula_1(h, otrezok[i][1], otrezok[i][2], priblijenie_2[2], otrezok[i][3], priblijenie_2[3]);
			priblijenie_3[2] = func(priblijenie_3[0], priblijenie_3[1]);
			priblijenie_3[3] = func_proiz_2(priblijenie_3[0], priblijenie_3[1], priblijenie_3[2]);
			//
			if ((fabs(priblijenie_1[1] - priblijenie_3[1]) < eps) && (fabs(priblijenie_1[2] - priblijenie_3[2]) < eps) && (fabs(priblijenie_1[3] - priblijenie_3[3]) < eps))
			{
				otrezok[i + 1][0] = priblijenie_3[0];
				otrezok[i + 1][1] = priblijenie_3[1];
				otrezok[i + 1][2] = priblijenie_3[2];
				otrezok[i + 1][3] = priblijenie_3[3];
				break;
			}
			else
			{
				h = h / 2;
				for (int i = 1; i < 3; i++)
				{
					otrezok[i][0] = otrezok[i - 1][0] + h;
					otrezok[i][1] = otrezok[i - 1][1] + delt(otrezok[i - 1][0], otrezok[i - 1][1], h);
					otrezok[i][2] = func(otrezok[i][0], otrezok[i][1]);
					otrezok[i][3] = func_proiz_2(otrezok[i][0], otrezok[i][1], otrezok[i][2]);
				}
			}
		}
	}
	return otrezok;
}
// дл€ 3 степени
double formula_1_3(double yim2, double yi, double yim1, double h, double y_3_i, double y_3_im1)
{
	return yim2 + 3.0 * (yi - yim1) + (h*h*h / 2) * (y_3_i + y_3_im1);
}
double formula_2_3(double yi, double y_1_ip1, double y_1_i, double y_2_ip1, double y_2_i, double y_3_ip1, double y_3_i,double h)
{
	return yi + h / 2 * (y_1_ip1 + y_1_i) - h*h / 10 * (y_2_ip1 - y_2_i) + h*h*h / 120 * (y_3_ip1 + y_3_i);
}
decltype(auto) method_2(double x0, double y0, double h, double eps, int n)
{
	vector<vector<double>> otrezok(n);
	for (auto& x : otrezok)
		x.resize(5);
	otrezok[0][0] = x0;
	otrezok[0][1] = y0;
	otrezok[0][2] = func(otrezok[0][0], otrezok[0][1]);
	otrezok[0][3] = func_proiz_2(otrezok[0][0], otrezok[0][1], otrezok[0][2]);
	otrezok[0][4] = func_proiz_3(otrezok[0][0], otrezok[0][1], otrezok[0][2], otrezok[0][3]);
	for (int i = 1; i < 3; i++)
	{
		otrezok[i][0] = otrezok[i - 1][0] + h;
		otrezok[i][1] = otrezok[i - 1][1] + delt(otrezok[i - 1][0], otrezok[i - 1][1], h);
		otrezok[i][2] = func(otrezok[i][0], otrezok[i][1]);
		otrezok[i][3] = func_proiz_2(otrezok[i][0], otrezok[i][1], otrezok[i][2]);
		otrezok[i][4] = func_proiz_3(otrezok[i][0], otrezok[i][1], otrezok[i][2], otrezok[i][3]);
	}
	vector<double>priblijenie_1(5);
	vector<double>priblijenie_2(5);
	vector<double>priblijenie_3(5);
	for (int i = 2; i < n - 1; i++)
	{
		while (true)
		{
			//
			priblijenie_1[0] = otrezok[i][0] + h;
			priblijenie_1[1] = formula_1_3(otrezok[i - 2][1], otrezok[i][1], otrezok[i - 1][1], h, otrezok[i][4], otrezok[i - 1][4]);
			priblijenie_1[2] = func(priblijenie_1[0], priblijenie_1[1]);
			priblijenie_1[3] = func_proiz_2(priblijenie_1[0], priblijenie_1[1], priblijenie_1[2]);
			priblijenie_1[4] = func_proiz_3(priblijenie_1[0], priblijenie_1[1], priblijenie_1[2], priblijenie_1[3]);
			//
			priblijenie_2[0] = otrezok[i][0] + h;
			priblijenie_2[1] = formula_2_3(otrezok[i][1], priblijenie_1[2], otrezok[i][2], priblijenie_1[3], otrezok[i][3], priblijenie_1[4], otrezok[i][4], h);
			priblijenie_2[2] = func(priblijenie_2[0], priblijenie_2[1]);
			priblijenie_2[3] = func_proiz_2(priblijenie_2[0], priblijenie_2[1], priblijenie_2[2]);
			priblijenie_2[4] = func_proiz_3(priblijenie_2[0], priblijenie_2[1], priblijenie_2[2], priblijenie_2[3]);
		
			//
			priblijenie_3[0] = otrezok[i][0] + h;
			priblijenie_3[1]= formula_2_3(otrezok[i][1], priblijenie_2[2], otrezok[i][2], priblijenie_2[3], otrezok[i][3], priblijenie_2[4], otrezok[i][4], h);
			priblijenie_3[2] = func(priblijenie_3[0], priblijenie_3[1]);
			priblijenie_3[3] = func_proiz_2(priblijenie_3[0], priblijenie_3[1], priblijenie_3[2]);
			priblijenie_3[4] = func_proiz_3(priblijenie_3[0], priblijenie_3[1], priblijenie_3[2], priblijenie_3[3]);
			if ((fabs(priblijenie_1[1] - priblijenie_3[1]) < eps) && (fabs(priblijenie_1[2] - priblijenie_3[2]) < eps) && (fabs(priblijenie_1[3] - priblijenie_3[3]) < eps) && (fabs(priblijenie_1[4] - priblijenie_3[4]) < eps))
			{
				otrezok[i + 1][0] = priblijenie_3[0];
				otrezok[i + 1][1] = priblijenie_3[1];
				otrezok[i + 1][2] = priblijenie_3[2];
				otrezok[i + 1][3] = priblijenie_3[3];
				otrezok[i + 1][4] = priblijenie_3[4];
				break;
			}
			else
			{
				h = h / 2;
				
				for (int i = 1; i < 3; i++)
				{
					otrezok[i][0] = otrezok[i - 1][0] + h;
					otrezok[i][1] = otrezok[i - 1][1] + delt(otrezok[i - 1][0], otrezok[i - 1][1], h);
					otrezok[i][2] = func(otrezok[i][0], otrezok[i][1]);
					otrezok[i][3] = func_proiz_2(otrezok[i][0], otrezok[i][1], otrezok[i][2]);
					otrezok[i][4]= func_proiz_3(otrezok[i][0], otrezok[i][1], otrezok[i][2], otrezok[i][3]);
				}
			}
		}
	}

	return otrezok;

}
int main()
{
	setlocale(LC_ALL, "Rus");
	double b;
	cout << "Ќачало отрезка: 1" << endl;
	cout << " онец отрезка: ";
	cin >> b;
	cout << "n=";
	int n;
	cin >> n;
	double a = 1;
	double h = (b - a) / n;
	double eps = 0.01;
	double y0 = 2;
	// метод высших пор€дков дл€ производной 2 степени
	vector<vector<double>>otrezok=method_2(a,y0,h,eps,n+1);
	// метод высших пор€дков дл€ производной 3 степени
	vector<vector<double>> gg = method_1(a, y0, h, eps, n+1);
	cout << endl;
	cout << "ћетод производных второго пор€дка" << endl;
	cout << "____________________________________________________________________________________________________________" << endl;
	cout << "|        x        |          y      |        y'       |         y''     |       f(x)     |     |y-f(x)|    |" << endl;
	cout << "|-----------------|-----------------|-----------------|-----------------|----------------|-----------------|" << endl;
	auto old=cout.setf(ios_base::fixed, ios_base::floatfield | ios_base::adjustfield);
	cout.precision(12);
	for (int i = 0; i < n; i++)
	{
		cout << "|";
		for (int j = 0; j < 4; j++)
		{
			cout.width(17);
			cout << gg[i][j] <<"|";
		}
		cout.width(16);
		cout << original_func(gg[i][0])<<"|";
		cout.width(17);
		cout<< fabs(original_func(gg[i][0]) - gg[i][1]) << "|";
		cout << endl;
	}
	cout << "|-----------------|-----------------|-----------------|-----------------|----------------|-----------------|" << endl;
	cout << endl;
	cout.setf(old);
	cout << "ћетод производных третьего пор€дка" << endl;
	cout << "______________________________________________________________________________________________________________________________" << endl;
	cout << "|        x        |          y      |        y'       |         y''     |        y'''     |       f(x)     |     |y-f(x)|    |" << endl;
	cout << "|-----------------|-----------------|-----------------|-----------------|-----------------|----------------|-----------------|" << endl;
	cout.setf(ios_base::fixed, ios_base::floatfield | ios_base::adjustfield);
	cout.precision(12);
	for (int i = 0; i < n; i++)
	{
		cout << "|";
		for (int j = 0; j < 5; j++)
		{
			cout.width(17);
			cout << otrezok[i][j] << "|";
		}
		cout.width(16);
		cout << original_func(otrezok[i][0]) << "|";
		cout.width(17);
		cout << fabs(original_func(otrezok[i][0]) - otrezok[i][1])<<"|";
		cout << endl;
	}
	system("pause");

}
