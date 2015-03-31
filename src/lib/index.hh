#ifndef  DUNE_INDEX_HH
# define DUNE_INDEX_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <iostream>
#include <fstream>

#include "dune/grid/common/mcmgmapper.hh" 
#include "dune/common/dynvector.hh" 

template<int dim> 
class Index 
{
	public:
		enum { dimension = dim };

		Index (const Dune::FieldVector<int, dim> & N) : N(N) { x = 0;};

		const Dune::FieldVector<int, dim> index () const {return x;};
		int index (int i) const {return x[i];};

		bool increaseIndex ()
		{
			int idim = dim - 1;
			while (true)
			{
				x[idim]++;
				if (x[idim] >= N[idim])
					idim--;
				else 
					break;
				if (idim < 0)
					return false;
			}

	//		std::cerr << "increaseIndex(): " << x << std::endl;

			for (int jdim = idim + 1; jdim < dim; jdim++)
				x[jdim] = 0;

			return true;
		}

	private:
		const Dune::FieldVector<int, dim> & N;
		Dune::FieldVector<int, dim> x;

};

void add (int a, int * res)
{
	*res = 10 + a;
}

#endif
