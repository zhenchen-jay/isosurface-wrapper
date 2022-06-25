#ifndef QEFQR_H
#define QEFQR_H

#include <stdio.h>
#include <vector>
using namespace std;

template <class Type, int n>
struct ArrayWrapper
{
	Type data [ n ] [ n ];
};


template <class Type, int n>
class QEFQR
{
public:
	Type data [ ( n + 1 ) * ( n + 2 ) / 2 ];

	void zero ( void )
	{
		for (int i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
			data [ i ] = 0;
	}

	void combineSelf ( Type *eqn );
	void combineSelf ( QEFQR<Type, n> &qef );

	void calcPoint (Type rvalue[], Type *w, ArrayWrapper<Type, n> &u ) const;
	void calcPoint (ArrayWrapper<Type, n> &a, Type *b, Type rvalue[], Type *w, ArrayWrapper<Type, n> &u ) const;
	void calcPointLower (ArrayWrapper<Type, n> &a, Type *b, Type rvalue[], Type *w, ArrayWrapper<Type, n-1> &u ) const;

	void build_ab (ArrayWrapper<Type, n> &a, Type *b) const;
};
















//** implementation


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

// this is the worst kind of hack, but the VC++ compiler can't do templates correctly
// In the case below, the compiler seems to be treating the variable "n" as a dynamic
// variable and not using the code transformation the template defines.
/*
template <class Type, int n>
void jacobi ( Type u[n][n], Type *d, Type v[n][n]);  <- Doesn't compile in VC++
*/

template <class Type, int n>
void jacobi ( ArrayWrapper<Type, n> &u, Type *d, ArrayWrapper<Type, n> &v )
{
	int j, iq, ip, i, k;
	Type tresh, theta, tau, t, sm, s, h, g, c;
	Type b [ n ];
	Type z [ n ];
	Type a [ n ] [ n ];

	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			a [ i ] [ j ] = u.data [ i ] [ j ];
		}
	}

	for ( ip = 0; ip < n; ip++ ) 
	{
		for ( iq = 0; iq < n; iq++ )
		{
			v.data [ ip ] [ iq ] = 0.0f;
		}
		v.data [ ip ] [ ip ] = 1.0f;
	}

	for ( ip = 0; ip < n; ip++ )
	{
		b [ ip ] = a [ ip ] [ ip ];
		d [ ip ] = b [ ip ];
		z [ ip ] = 0.0f;
	}

	for ( i = 1; i <= 50; i++ )
	{
		sm = 0.0f;
		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ )
			{
				sm += (Type)fabs ( a [ ip ] [ iq ] );
			}
		}

		if ( sm == 0.0f )
		{
			// sort the stupid things and transpose

			// transpose
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					a [ i ] [ j ] = v.data [ j ] [ i ];
				}
			}

			// stupid sort n^2... however n should always be small
			// bubble sort
			// largest first
			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n - i - 1; j++ )
				{
					if ( fabs ( d [ j ] ) < fabs ( d [ j + 1 ] ) )
					{
						sm = d [ j ];
						d [ j ] = d [ j + 1 ];
						d [ j + 1 ] = sm;

						for ( k = 0; k < n; k++ )
						{
							sm = a [ j ] [ k ];
							a [ j ] [ k ] = a [ j + 1 ] [ k ];
							a [ j + 1 ] [ k ] = sm;
						}
					}
				}
			}

			for ( i = 0; i < n; i++ )
			{
				for ( j = 0; j < n; j++ )
				{
					v.data [ i ] [ j ] = a [ i ] [ j ];
				}
			}

			return;
		}

		if ( i < 4 )
		{
			tresh = 0.2f * sm / ( n * n );
		}
		else
		{
			tresh = 0.0f;
		}

		for ( ip = 0; ip < n - 1; ip++ )
		{
			for ( iq = ip + 1; iq < n; iq++ ) 
			{
				g = 100.0f * (Type)fabs ( a [ ip ] [ iq ] );
				if ( i > 4 && ( fabs ( d [ ip ] ) + g ) == fabs ( d [ ip ] )
					&& ( fabs ( d [ iq ] ) + g ) == fabs ( d [ iq ] ) )
				{
					a [ ip ] [ iq ] = 0.0f;
				}
				else
				{
					if ( fabs ( a [ ip ] [ iq ] ) > tresh )
					{
						h = d [ iq ] - d [ ip ];
						if ( ( fabs ( h ) + g ) == fabs ( h ) )
						{
							t = ( a [ ip ] [ iq ] ) / h;
						}
						else
						{
							theta = 0.5f * h / ( a [ ip ] [ iq ] );
							t = 1.0f / (Type)( fabs ( theta ) + sqrt ( 1.0f + theta * theta ) );
							if ( theta < 0.0f ) 
							{
								t = -1.0f * t;
							}
						}

						c = 1.0f / (Type)sqrt ( 1 + t * t );
						s = t * c;
						tau = s / ( 1.0f + c );
						h = t * a [ ip ] [ iq ];
						z [ ip ] -= h;
						z [ iq ] += h;
						d [ ip ] -= h;
						d [ iq ] += h;
						a [ ip ] [ iq ] = 0.0f;
						for ( j = 0; j <= ip - 1; j++ )
						{
							ROTATE ( a, j, ip, j, iq )
						}
						for ( j = ip + 1; j <= iq - 1; j++ )
						{
							ROTATE ( a, ip, j, j, iq )
						}
						for ( j = iq + 1; j < n; j++ )
						{
							ROTATE ( a, ip, j, iq, j )
						}
						for ( j = 0; j < n; j++ )
						{
							ROTATE ( v.data, j, ip, j, iq )
						}
					}
				}
			}
		}

		for ( ip = 0; ip < n; ip++ )
		{
			b [ ip ] += z [ ip ];
			d [ ip ] = b [ ip ];
			z [ ip ] = 0.0f;
		}
	}

	printf ( "too many iterations in jacobi\n" );
	exit ( 1 );
}

template <class Type, int n>
void matInverse ( ArrayWrapper<Type, n> &mat, ArrayWrapper<Type, n> &rvalue, Type tolerance, Type *w, ArrayWrapper<Type, n> &u )
{
	// there is an implicit assumption that mat is symmetric and real
	// U and V in the SVD will then be the same matrix whose rows are the eigenvectors of mat
	// W will just be the eigenvalues of mat
	int i, j, k;

	jacobi<Type, n> ( mat, w, u );

	// pseudo invert eigenvalues
	for ( i = 1; i < n; i++ )
	{
		if ( fabs ( w [ i ] / w [ 0 ] ) < tolerance )
			w [ i ] = 0;
		else
			w [ i ] = 1.0f / w [ i ];
	}
	w [ 0 ] = 1.0f / w [ 0 ];

	// multiply back to get inverted matrix
	for ( i = 0; i < n; i++ )
	{
		for ( j = 0; j < n; j++ )
		{
			rvalue.data [ i ] [ j ] = 0;
			for ( k = 0; k < n; k++ )
			{
				rvalue.data [ i ] [ j ] += w [ k ] * u.data [ k ] [ i ] * u.data [ k ] [ j ];
			}
		}
	}
}






















template <class Type, int n>
void QEFQR<Type, n>::combineSelf ( Type *eqn )
{
	int i, j, index;
	Type a, b, mag, temp;
	int rvalueN = 0;
	Type scale = 0;

	eqn [ n ] *= -1;

	for ( i = 0; i < n + 1; i++ )
	{
		index = ( ( 2 * n + 3 - i ) * i ) / 2;

		a = data [ index ];
		b = eqn [ i ];

		if ( fabs ( a ) > 0 || fabs ( b ) > 0 )
		{
			rvalueN = i + 1;

			mag = (Type)sqrt ( a * a + b * b );
			a /= mag;
			b /= mag;

			for ( j = 0; j < n + 1 - i; j++ )
			{
				temp = a * data [ index + j ] + b * eqn [ j + i ];
				eqn [ j + i ] = b * data [ index + j ] - a * eqn [ j + i ];
				data [ index + j ] = temp;
			}
		}
	}

	eqn [ n ] *= -1;

	// print the equation
	//fprintf(trace_file, "qef data = ");
	//for ( i = 0; i < ( n + 1 ) * ( n + 2 ) / 2; i++ )
	//{
	//	fprintf(trace_file, "%I64x ", *(__int64*)&data [ i ]);
	//	//fprintf(trace_file, "%g ", data [ i ]);
	//}
	//fprintf(trace_file, "\n");
}

template <class Type, int n>
void QEFQR<Type, n>::build_ab (ArrayWrapper<Type, n> &a, Type *b) const
{
	int i, j, k;

	for ( i = 0; i < n; i++ )
	{
		b [ i ] = 0;
		for ( j = 0; j < n; j++ )
			a.data [ i ] [ j ] = 0;
	}

	for ( k = 0; k < n; k++ )
	{
		const int index = ( ( 2 * n + 3 - k ) * k ) / 2;
		for ( i = k; i < n; i++ )
		{
			for ( j = k; j < n; j++ )
				a.data [ i ] [ j ] += data [ index + i - k ] * data [ index + j - k ];
		}
	}

	for ( k = 0; k < n; k++ )
	{
		for ( i = 0; i <= k; i++ )
		{
			const int index = ( ( 2 * n + 3 - i ) * i ) / 2;
			b [ k ] += data [ index + k - i ] * data [ index + n - i ];
		}
	}
}

template <class Type, int n>
void QEFQR<Type, n>::calcPoint (Type rvalue[], Type *w, ArrayWrapper<Type, n> &u ) const
{
	ArrayWrapper<Type, n> a;
	Type b [ n ];
	build_ab(a, b);

	calcPoint(a,b, rvalue, w,u);
}

template <class Type, int n>
void QEFQR<Type, n>::calcPoint (ArrayWrapper<Type, n> &a, Type *b, Type rvalue[], Type *w, ArrayWrapper<Type, n> &u ) const
{
	ArrayWrapper<Type, n> inv;
	::matInverse<Type, n> (a, inv, 1e-6, w, u);

	for (int i = 0; i < n; i++)
	{
		rvalue [ i ] = 0;
		for (int j = 0; j < n; j++)
			rvalue[i] += inv.data[j][i]*b[j];
	}
}


template <class Type, int n>
void QEFQR<Type, n>::calcPointLower (ArrayWrapper<Type, n> &a, Type *b, Type rvalue[], Type *w, ArrayWrapper<Type, n-1> &u) const
{
	ArrayWrapper<Type, n-1> aa;

	for (int i = 0; i < n-1; i++)
		for (int j = 0; j < n-1; j++)
			aa.data[i][j] = a.data[i][j];

	ArrayWrapper<Type, n-1> inv;
	::matInverse<Type, n-1> (aa, inv, 1e-6, w, u);

	for (int i = 0; i < n-1; i++)
	{
		rvalue [ i ] = 0;
		for (int j = 0; j < n-1; j++)
			rvalue[i] += inv.data[j][i]*b[j];
	}
}



#endif