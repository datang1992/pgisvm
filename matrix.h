

#ifndef __MATRIX__H___
#define __MATRIX__H___


class matrix
{	
	public:
		int m, n;
		double **num;
	
		matrix();
		matrix(int, int); //construct an identity matrix if the second index is non-positive, the first index will be the dimension of the identity matrix.
		matrix(const matrix&);

		matrix& operator= (const matrix&);
		matrix operator+ (const matrix&) const;
		matrix operator- () const;
		matrix operator- (const matrix&) const;
		matrix operator* (const matrix&) const;
		matrix operator* (const double&) const;
		matrix inverse () const;
		matrix transpose () const;
		matrix cholesky () const;
		double det() const;
        double trace() const;

		~matrix();
};







#endif
