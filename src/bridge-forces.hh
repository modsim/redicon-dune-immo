/*  testBridgeForce.cc  2015-03-04 test Bridge class
 *
 * Copyright (C) 2015 Svyatoslav Kondrat (Valiska)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <iostream>
#include "bridge.hh"

#define pow2(x) 	 ((x)*(x))

template <typename ct, int dim, typename ft, int fdim>
class ForceElectrostatic
{
	public: 
		enum {fdimension = fdim * dim};
		typedef double ftype;

		ForceElectrostatic (const Dune::FieldVector<ft, fdim> & A, ct d) : A(A), d(d) {};
		Dune::FieldVector<ftype, fdimension> 
			operator () (const Dune::FieldVector<ct, dim> & x, const Dune::FieldVector<ct, dim> & y) const 
			{
				ct r3 = 0.;
				ct dx[dim] = {0,0,0};
				for (int i = 0; i < dim - 1; i++)
				{
					dx[i] = x[i] - y[i];
					r3 += pow2 (dx[i]);
				}
				// last coordinate (plate) is molecular diameter away
				dx[dim-1] = x[dim-1] - y[dim-1] + d;
				r3 += pow2 (dx[dim-1]); // r^2

				r3 *= sqrt (r3);      // r^3
				Dune::FieldVector<ftype, fdimension> val (1./r3);
				for (int a = 0; a < fdim; a++)
					for (int i = 0; i < dim; i++)
					{
						val[a*dim + i] *= - A[a] * dx[i]; 
					}
				return val;
			};
	private:
		const Dune::FieldVector<ft, fdim> & A;
		const ct d;
};

template <typename ct, int dim, typename ft, int fdim>
class ForceLJ
{
	public: 
		enum {fdimension = fdim * dim};
		typedef double ftype;

		ForceLJ (const Dune::FieldVector<ft, fdim> & A, ct d) : A(A), d(d) {};
		Dune::FieldVector<ftype, fdimension> 
			operator () (const Dune::FieldVector<ct, dim> & x, const Dune::FieldVector<ct, dim> & y) const 
			{
				ct r2 = 0.;
				ct dx[dim] = {0,0,0};
				for (int i = 0; i < dim - 1; i++)
				{
					dx[i] = x[i] - y[i];
					r2 += pow2 (dx[i]); // r^2
				}
				// last coordinate (plate) is molecular diameter away
				dx[dim-1] = x[dim-1] - y[dim-1] + d;
				r2 += pow2 (dx[dim-1]); // r^2

				ct r4 = pow2(r2); // r^4
				ct r8 = pow2(r4); // r^8
				ct r8_1 = 1./r8;
				//if ( fabs(r8_1) < 10.E-8 )
				//	r8_1 = 0.;
					
				Dune::FieldVector<ftype, fdimension> val (r8_1);
				for (int a = 0; a < fdim; a++)
					for (int i = 0; i < dim; i++)
					{
						val[a*dim + i] *= - A[a] * dx[i]; 
					}
				return val;
			};
	private:
		const Dune::FieldVector<ft, fdim> & A;
		const ct d;
};


