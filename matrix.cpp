

#include "matrix.h"
#include <math.h>
#include <iostream>


using namespace std;


matrix::matrix()
{
	m = 0;
	n = 0;
	num = NULL;
}

matrix::matrix(int M, int N)
{	
	m = M;
	n = N;

	num = new double*[m];

	if(N > 0)
	{
		for(int i = 0; i < m; i++)
		{
			num[i] = new double[n];
			for(int j = 0; j < n; j++)
				num[i][j] = 0;
		}
	}
	else
	{
		if(N > -100)
		{
			n = m;
			for(int i = 0; i < m; i++)
			{
				num[i] = new double[n];
				for(int j = 0; j < n; j++)
					num[i][j] = 0;
				num[i][i]  = 1;
			}
		}
		else if(N > -10000)
		{
			n = 1;
			for(int i = 0; i < m; i++)
			{
				num[i] = new double[1];
				num[i][0] = 0;
			}
		}
	}
}


matrix::matrix (const matrix& mat)
{

	m = mat.m;
	n = mat.n;
	

	if(mat.num)
	{
		num = new double*[m];
		for(int i = 0; i < m; i++)
			if(mat.num[i])
			{
				num[i] = new double[n];
				for(int j = 0; j < n; j++)
					num[i][j] = mat.num[i][j];
			}
	}
}

matrix& matrix::operator= (const matrix& mat)
{
	m = mat.m;
	n = mat.n;
	
	if(num)
	{
		
        for(int i = 0; i < m; i++)
			if(num[i])
				delete[] num[i];
		delete[] num;
        
        
	}

	if(mat.num)
	{
		num = new double*[m];
		for(int i = 0; i < m; i++)
			if(mat.num[i])
			{
				num[i] = new double[n];
				for(int j = 0; j < n; j++)
					num[i][j] = mat.num[i][j];
			}
	}
	return (*this);
}


matrix matrix::operator+ (const matrix& mat) const
{
	if(m != mat.m || n != mat.n)
	{
		cout << "Error: Add two matrices with different size! A: " << m << " * " << n << ", B: " << mat.m << " * " << n <<  "!" << endl;
		return matrix();
	}
	
	matrix ret = matrix(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[i][j] = num[i][j] + mat.num[i][j];

	return ret;
}

matrix matrix::operator- () const
{
	matrix  ret = matrix(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[i][j] = -num[i][j];
	return ret;
}


matrix matrix::operator- (const matrix& mat) const
{
	if(m != mat.m || n != mat.n)
	{
		cout << "Error: Sub matrix B from A,  with different size! A: " << m << " * " << n << ", B: " << mat.m << " * " << n <<  "!" << endl;
		return matrix();
	}

	matrix ret = matrix(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[i][j] = num[i][j] - mat.num[i][j];

	return ret;
}

matrix matrix::transpose () const
{
	matrix  ret = matrix(n, m);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[j][i] = num[i][j];
	return ret;
}

matrix matrix::operator* (const matrix& mat) const
{
	if(n != mat.m)
	{
		cout << "Error: Calculate A * B, with A.col = " << n << " while B.row = " << mat.m << "!" << endl;
		return matrix();
	}

	matrix ret = matrix(m, mat.n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < mat.n; j++)
		{
			ret.num[i][j] = 0;
			for(int k = 0; k < n; k++)
				ret.num[i][j] += num[i][k] * mat.num[k][j];
		}
	return ret;
}


matrix matrix::operator* (const double& number) const
{
	matrix ret = matrix(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[i][j] = number * num[i][j];
	return ret;
}


matrix matrix::inverse() const
{
	if(m != n)
	{
		cout << "Error: Calculate inverse of a matrix with is not a square matrix!" << endl;
		return matrix();
	}

	double **tmp = new double*[n];

	for(int i = 0; i < n; i++)
	{
		tmp[i] = new double[n * 2];
		for(int j = 0; j < n; j++)
			tmp[i][j] = num[i][j];
		for(int j = n; j < n * 2; j++)
			tmp[i][j] = 0;
		tmp[i][i + n] = 1;
	}
	
	for(int i = 0; i < n; i++)
	{
		double maxv = fabs(tmp[i][i]);
		int maxn = i;
		for(int j = i + 1; j < n; j++)
			if(fabs(tmp[j][i]) > maxv)
			{
				maxv = fabs(tmp[j][i]);
				maxn = j;
			}
		
		if(maxv < 1e-6)
		{
			cout << "Error: Try to inverse an un-inversible matrix!" << endl;
			return matrix();
		}

		for(int j = i; j < n * 2; j++)
		{
			double t = tmp[i][j];
			tmp[i][j] = tmp[maxn][j];
			tmp[maxn][j] = t;
		}
		
		for(int j = n * 2 - 1; j >= i; j--)
			tmp[i][j] /= tmp[i][i];
		
		
		for(int j = i + 1; j < n; j++)
			for(int k = n * 2 - 1; k >= i; k--)
				tmp[j][k] -= tmp[j][i] * tmp[i][k];

	}

	for(int i = n - 1; i >= 0; i--)
		for(int j = i - 1; j >= 0; j--)
			for(int k = n * 2 - 1; k >= i; k--)
				tmp[j][k] -= tmp[j][i] * tmp[i][k];

	matrix ret = matrix(m, n);
	for(int i = 0; i < m; i++)
		for(int j = 0; j < n; j++)
			ret.num[i][j] = tmp[i][j + n];

    
	for(int i = 0; i < m; i++)
		delete[] tmp[i];
	
	delete[] tmp;

    
	return ret;
}

double matrix::det() const
{
	if(m != n)
	{
		cout << "Error: Calculate determinant of a matrix with is not a square matrix!" << endl;
		return 0;
	}

	matrix tmp = matrix(*this);
    double ret = 1;
	for(int i = 0; i < n; i++)
	{
		double maxv = 0;
        int maxn = i;
		for(int j = i; j < n; j++)
			if(fabs(tmp.num[j][i]) > maxv)
			{
				maxv = fabs(tmp.num[j][i]);
				maxn = j;
			}
		if(maxv < 1e-6)
			return 0;
		if(maxn != i)
		  ret *= -1;
		for(int j = i; j < n; j++)
			swap(tmp.num[i][j], tmp.num[maxn][j]);
		for(int j = i + 1; j < n; j++)
		{
			double ratio = -tmp.num[j][i] / tmp.num[i][i];
			for(int k = i; k < n; k++)
				tmp.num[j][k] += ratio * tmp.num[i][k];
		}
		ret *= tmp.num[i][i];
	}
	return ret;
}

double matrix::trace() const
{
	if(m != n)
	{
		cout << "Error: Calculate trace of a matrix with is not a square matrix!" << endl;
		return 0;
	}
    
    double ret = 0;
    for (int i = 0; i < m; i++) {
        ret += num[i][i];
    }
	return ret;
}



matrix matrix::cholesky() const
{
	if(m != n)
	{
		cout << "Error: Try to find the Cholosky decomposition of a non-square matrix!" << endl;
		return matrix();
	}

	matrix ret = matrix(m, n);
	
	for(int i = 0; i < n; i++)
	{
		for(int j = 0; j <= i; j++)
			ret.num[i][j] = num[i][j];
		for(int j = i + 1; j < n; j++)
			ret.num[i][j] = 0;
	}

	for(int i = 0; i < n; i++)
	{
		for(int k = 0; k < i; k++)
			ret.num[i][i] -= ret.num[i][k] * ret.num[i][k];
		if(num[i][i] < 0)
		{
			cout << "Error: Try to find the Cholosky decompostion of a non-positive-semi-definate matrix!" << endl;
			return matrix();
		}
		ret.num[i][i] = sqrt(ret.num[i][i]);
		
		for(int j = i + 1; j < n; j++)
		{
			for(int k = 0; k < i; k++)
				ret.num[j][i] -= ret.num[j][k] * ret.num[i][k];
			ret.num[j][i] /= ret.num[i][i];
		}
	}

	return ret;
}


matrix::~matrix()
{
	
    if(num)
	{
		for(int i = 0; i < m; i++)
			delete[] num[i];
		delete[] num;
	}
    

}


