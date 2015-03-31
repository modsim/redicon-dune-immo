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

template <typename ct, int dim>
class Force
{
	public: 
		enum {fdimension = 3};
		typedef double ftype;

		Force (double a) : a(a) {};
		Dune::FieldVector<ftype, fdimension> 
			operator () (const Dune::FieldVector<ct, dim> & x, const Dune::FieldVector<ct, dim> & y) const 
			{
				double r2 = 0;
				double dx[3] = {0,0,0};
				for (int i = 0; i < dim - 1; i++)
				{
					dx[i] = x[i] - y[i];
					r2 += pow2 (dx[i]);
				}
				dx[dim-1] = x[dim-1] - y[dim-1] + 2 * a;
				r2 += pow2 (dx[dim-1]);

				r2 *= sqrt (r2);
				Dune::FieldVector<ftype, fdimension> val (1/r2);
				for (int i = 0; i < dim; i++)
				{
					val[i] *= - dx[i]; 
				}
				return val;
			};
	private:
		const double a;
};

#include "../bridge-forces.hh"

int main (int argc, char ** argv) 
{
	Dune::FieldVector<double, 3> x0;
	x0[0] = 15.; x0[1] = 15.; x0[2] = 0.;

	Dune::FieldVector<double, 2> H;
	H[0] = 10.; H[1] = 10.;

	int p = 40;

	typedef Force<double, 3>  Function;
	const Function & F(1.);
	RectBridge<double, 3, Function> B = RectBridge<double, 3, Function> (x0, H, F, 50);

	Dune::FieldVector<double, 1> Q (-.1);

	typedef ForceElectrostatic<double, 3, double, 1>  FunctionEl;
	FunctionEl Fel(Q, 0.5);
	RectBridge<double, 3, FunctionEl> Bel = RectBridge<double, 3, FunctionEl> (x0, H, Fel, p);

	typedef ForceLJ<double, 3, double, 1>  FunctionLJ6;
	FunctionLJ6 Flj6(Q, 0.5);
	RectBridge<double, 3, FunctionLJ6> Blj6 = RectBridge<double, 3, FunctionLJ6> (x0, H, Flj6, p);

	Dune::FieldVector<double, 3>  x (15.);
	//x[0] = 5.;
	x[2] = 0.0;
	double step = 0.2;

	for (int i = 0; i < 100; i++)
	{
		std::cout << x << "  " << Bel(x) << "   " << Blj6(x) << std::endl;
		x[2] += step;
	}

	return 1;
}

