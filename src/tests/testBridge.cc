/*  testBridge.cc  2015-03-03 test BridgePotential class
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
//#include <random>

#include "bridge.hh"

template <typename ct, int dim>
class Integrand
{
	public: 
		enum {fdimension = 1};
		typedef double ftype;

		Integrand (int i) {};
		Dune::FieldVector<ftype, fdimension> 
			operator () (const Dune::FieldVector<ct, dim> & x, const Dune::FieldVector<ct, dim> & y) const 
			{
				double r2 = 0;
				for (int i = 0; i < dim; i++)
				{
					double dx = x[i] - y[i];
					r2 += pow2 (dx);
				}
				Dune::FieldVector<ftype, fdimension> val (1/sqrt(r2));
				return val;
			};
};

int main (int argc, char ** argv) 
{
	Dune::FieldVector<double, 3> x0;
	x0[0] = 15.; x0[1] = 15.; x0[2] = 0.;

	Dune::FieldVector<double, 2> H;
	H[0] = 10.; H[1] = 10.;

	typedef Integrand<double, 3>  Function;
	const Function & F(0);
	RectBridge<double, 3, Function> B = RectBridge<double, 3, Function> (x0, H, F, 50);

	Dune::FieldVector<double, 3>  x (10.);
	x[2] = 0.0;
	double step = 0.2;

	for (int i = 0; i < 100; i++)
	{
		std::cout << x << "  " << B(x) << std::endl;
		x[2] += step;
	}

	return 1;
}

