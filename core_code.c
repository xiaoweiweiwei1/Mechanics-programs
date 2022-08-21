#include <stdio.h>
#define v 0.2
#define E 1

void Kmatrix(double xi, double yi, double xj, double yj, double xk, double yk, double K[][6])
{
	double bi = yj - yk, bj = yk - yi, bk = yi-yj, ci = -xj + xk, cj = -xk + xi, ck = -xi + xj;
	double b[3] = { bi,bj,bk }, c[3] = { ci,cj,ck }, k[6][6], A = ((yk - yi) * (xj - xi) - (yi - yj) * (xi - xk)) / 2;
	int row, col;
	for (row = 0; row < 5; row += 2)
	{
		for (col = 0; col < 5; col += 2)
			k[row][col] = b[row / 2] * b[col / 2] + (1 - 0.2) / 2 * c[row / 2] * c[col / 2];
		for (col = 1; col < 6; col += 2)
			k[row][col] = v * b[row / 2] * c[(col - 1) / 2] + (1 - v) / 2 * c[row / 2] * b[(col - 1) / 2];
	}
	for (row = 1; row < 6; row += 2)
	{
		for (col = 0; col < 5; col += 2)
			k[row][col] = v * c[(row - 1) / 2] * b[col / 2] + (1 - v) / 2 * b[(row - 1) / 2] * c[col / 2];
		for (col = 1; col < 6; col += 2)
			k[row][col] = c[(row - 1) / 2] * c[(col - 1) / 2] + (1 - v) / 2 * b[(row - 1) / 2] * b[(col - 1) / 2];
	}
	for (row = 0; row < 6; row++)
		for (col = 0; col < 6; col++)
			K[row][col] = k[row][col] * E / (4 * (1 - v * v) * A);

}

void Bmatrix(double xi, double yi, double xj, double yj, double xk, double yk, double B[][6])
{
	double bi = yj - yk, bj = yk - yi, bk = yi - yj, ci = -xj + xk, cj = -xk + xi, ck = -xi + xj;
	double b[3][6] = { {bi,0,bj,0,bk,0},{0,ci,0,cj,0,ck},{ci,bi,cj,bj,ck,bk} }, A = ((yk - yi) * (xj - xi) - (yi - yj) * (xi - xk)) / 2;
	int row, col;
	for (row = 0; row < 3; row++)
		for (col = 0; col < 6; col++)
			B[row][col] = b[row][col] / (2 * A);
}

void Smatrix(double xi, double yi, double xj, double yj, double xk, double yk, double S[][6])
{
	double bi = yj - yk, bj = yk - yi, bk = yi - yj, ci = -xj + xk, cj = -xk + xi, ck = -xi + xj, A = ((yk - yi) * (xj - xi) - (yi - yj) * (xi - xk)) / 2;
	double s[3][6] = { {bi,v * ci,bj,v * cj,bk,v * ck},{v * bi,ci,v * bj,cj,v * bk,ck},{(1 - v) * ci / 2,(1 - v) * bi / 2,(1 - v) * cj / 2,(1 - v) * bj / 2,(1 - v) * ck / 2,(1 - v) * bk / 2} };
	int row, col;
	for (row = 0; row < 3; row++)
		for (col = 0; col < 6; col++)
			S[row][col] = s[row][col] * E / ((2 * A) * (1 - v * v));
}



void main() 
{
	double K1[6][6] = { 0 }, K2[6][6] = { 0 };
	double delta[8] = { 1,1,1,1,1,1,1,1 }, F[8] = { 0,0,0,0,100,0,0,0 }, a[8][8] = { 0 }, b[8] = { 0,0,0,0,100,0,0,0 };
	int i, j, k, n = 8;
	double sum = 0;
	double B1[3][6] = { 0 }, B2[3][6] = { 0 }, S1[3][6] = { 0 }, S2[3][6] = { 0 };
	double delta1[6] = { 0 }, delta2[6] = { 0 };
	double x1 = 0,y1 = 0,x2 = 0,y2 = 1,x3 = 1,y3 = 1,x4 = 1,y4 = 0;
	int d[8] = { 1,1,1,1,1,1,1,1 };
	double tau1, tau2, sigma1x, sigma1y, sigma2x, sigma2y, gamma1, gamma2, epsilon1x, epsilon1y,epsilon2x,epsilon2y;
	Kmatrix(x1, y1, x4, y4, x2, y2, K1);
	Kmatrix(x4, y4, x3, y3, x2, y2, K2);
	Bmatrix(x1, y1, x4, y4, x2, y2, B1);
	Bmatrix(x4, y4, x3, y3, x2, y2, B2);
	Smatrix(x1, y1, x4, y4, x2, y2, S1);
	Smatrix(x4, y4, x3, y3, x2, y2, S2);

	//a[8][8]
	for (i = 0; i < 2; i++)
		for (j = 0; j < 2; j++)
			a[i][j] = K1[i][j];
	for (i = 0; i < 2; i++)
		for (j = 2; j < 4; j++)
			a[i][j] = K1[i][j + 2];
	for (i = 0; i < 2; i++)
		for (j = 6; j < 8; j++)
			a[i][j] = K1[i][j - 4];
	for (i = 2; i < 4; i++)
		for (j = 0; j < 2; j++)
			a[i][j] = K1[i + 2][j];
	for (i = 2; i < 4; i++)
		for (j = 2; j < 4; j++)
			a[i][j] = K1[i + 2][j + 2] + K2[i + 2][j + 2];
	for (i = 2; i < 4; i++)
		for (j = 4; j < 6; j++)
			a[i][j] = K2[i + 2][j - 2];
	for (i = 2; i < 4; i++)
		for (j = 6; j < 8; j++)
			a[i][j] = K1[i + 2][j - 4] + K2[i + 2][j - 6];
	for (i = 4; i < 6; i++)
		for (j = 2; j < 4; j++)
			a[i][j] = K2[i - 2][j + 2];
	for (i = 4; i < 6; i++)
		for (j = 4; j < 6; j++)
			a[i][j] = K2[i - 2][j - 2];
	for (i = 4; i < 6; i++)
		for (j = 6; j < 8; j++)
			a[i][j] = K2[i - 2][j - 6];
	for (i = 6; i < 8; i++)
		for (j = 0; j < 2; j++)
			a[i][j] = K1[i - 4][j];
	for (i = 6; i < 8; i++)
		for (j = 2; j < 4; j++)
			a[i][j] = K1[i - 4][j + 2] + K2[i - 6][j + 2];
	for (i = 6; i < 8; i++)
		for (j = 4; j < 6; j++)
			a[i][j] = K2[i - 6][j - 2];
	for (i = 6; i < 8; i++)
		for (j = 6; j < 8; j++)
			a[i][j] = K1[i - 4][j - 4] + K2[i - 6][j - 6];

	//setting 0
	d[0] = 0; d[1] = 0; d[7] = 0;

	for (i = 0; i < 8; i++)
		if (d[i] == 0)
		{
			for (j = 0; j < 8; j++)
				a[i][j] = 0;
			for (j = 0; j < 8; j++)
				a[j][i] = 0;
			a[i][i] = 1;
		}

	//Jordan
	for (k = 0; k < n; k++)
	{
		for (j = k + 1; j < n; j++)
			a[k][j] = a[k][j] / a[k][k];
		b[k] = b[k] / a[k][k];
		for (i = 0; i < n; i++)
		{
			if (i == k) continue;
			for (j = k + 1; j < n; j++)
				a[i][j] = a[i][j] - a[i][k] * a[k][j];
				b[i] = b[i] - a[i][k] * b[k];
		}	
	}

	//epsilon
	delta1[0] = b[0]; delta1[1] = b[1]; delta1[2] = b[6]; delta1[3] = b[7]; delta1[4] = b[2]; delta1[5] = b[3];
	delta2[0] = b[6]; delta2[1] = b[7]; delta2[2] = b[4]; delta2[3] = b[5]; delta2[4] = b[2]; delta2[5] = b[3];

	sum = 0;
	for (i = 0, j = 0; j < 6; j++)
		sum = sum + B1[i][j] * delta1[j];
	epsilon1x = sum;

	sum = 0;
	for (i = 1, j = 0; j < 6; j++)
		sum = sum + B1[i][j] * delta1[j];
	epsilon1y = sum;

	sum = 0;
	for (i = 2, j = 0; j < 6; j++)
		sum = sum + B1[i][j] * delta1[j];
	gamma1 = sum;

	sum = 0;
	for (i = 0, j = 0; j < 6; j++)
		sum = sum + B2[i][j] * delta2[j];
	epsilon2x = sum;

	sum = 0;
	for (i = 1, j = 0; j < 6; j++)
		sum = sum + B2[i][j] * delta2[j];
	epsilon2y = sum;

	sum = 0;
	for (i = 2, j = 0; j < 6; j++)
		sum = sum + B2[i][j] * delta2[j];
	gamma2 = sum;

	//sigma
	sum = 0;
	for (i = 0, j = 0; j < 6; j++)
		sum = sum + S1[i][j] * delta1[j];
	sigma1x = sum;

	sum = 0;
	for (i = 1, j = 0; j < 6; j++)
		sum = sum + S1[i][j] * delta1[j];
	sigma1y = sum;

	sum = 0;
	for (i = 2, j = 0; j < 6; j++)
		sum = sum + S1[i][j] * delta1[j];
	tau1 = sum;

	sum = 0;
	for (i = 0, j = 0; j < 6; j++)
		sum = sum + S2[i][j] * delta2[j];
	sigma2x = sum;

	sum = 0;
	for (i = 1, j = 0; j < 6; j++)
		sum = sum + S2[i][j] * delta2[j];
	sigma2y = sum;

	sum = 0;
	for (i = 2, j = 0; j < 6; j++)
		sum = sum + S2[i][j] * delta2[j];
	tau2 = sum;

printf("%lf %lf %lf %lf %lf %lf", epsilon1x, epsilon1y,gamma1, epsilon2x, epsilon2y,gamma2);
printf("\n%lf %lf %lf %lf %lf %lf", sigma1x, sigma1y, tau1, sigma2x, sigma2y, tau2);
}