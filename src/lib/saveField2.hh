#ifndef  DUNE_SAVE_FIELD2_HH
# define DUNE_SAVE_FIELD2_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <iostream>
#include <fstream>

#include "dune/grid/common/mcmgmapper.hh" 
#include "dune/common/dynvector.hh" 
#include "index.hh"
#include "evalField.hh"

/*
template<int dim> bool increaseIndex (Dune::FieldVector<int, dim> * index,
		const Dune::FieldVector<int, dim> & N)
{
	int idim = dim-1;
	while (true)
	{
		(*index)[idim]++;
		if ((*index)[idim] >= N[idim])
			idim--;
		else 
			break;
		if (idim < 0)
			return false;
	}

	//std::cerr << "increaseIndex(): " << *index << std::endl;

	for (int jdim = idim + 1; jdim < dim; jdim++)
		(*index)[jdim] = 0;
	
	return true;
}
*/

template<class Grid, typename ft, int fdim> void saveField2 (const Grid & g, const Dune::DynamicVector<ft> & F, 
		const Dune::FieldVector<typename Grid::ctype, Grid::dimension> x0,
		const Dune::FieldVector<typename Grid::ctype, Grid::dimension> l,
		const Dune::FieldVector<int, Grid::dimension> N,
		std::ostream * stream)
{
	const int dim = Grid::dimension;
	typedef typename Grid::LeafGridView GridView;
	const GridView & gv = g.leafView();
		
	Dune::FieldVector<typename GridView::ctype, GridView::dimension> dx;
	for (int idim = 0; idim < dim; idim++)
		if (N[idim])
			dx[idim] = l[idim] / ((double) N[idim] - 1);
		else
			dx[idim] = 0.;
#ifdef DEBUG
	std::cerr << "dx=" << dx << std::endl;
#endif

	Dune::FieldVector<typename GridView::ctype, GridView::dimension> x = x0;
	Index<dim> ind (N);

	typedef evalField<GridView, double, fdim>  EvalField;
	EvalField f (gv, F);

	int nldim = dim - 1; // new line index
	for (int idim = dim - 1; idim >=0; idim--)
		if ( (N[idim] != 0) && (l[idim] != 0.0) )
		{
			nldim = idim; 
			break;
		}

	//std::cerr << "nldim = " << nldim << std::endl;

	for (int idim = 0; idim < dim; idim++)
		if ( (N[idim] == 0) || (l[idim] == 0.0) )
			*stream << "# Fixed coordinate: x"<< idim << "=" << x[idim] << std::endl;

	while (true)
	{
		for (int idim = 0; idim < dim; idim++)
			if ( (N[idim]!=0) && (l[idim] != 0.0) )
				*stream << x[idim] << "   ";

		*stream << f(x) << std::endl;

		if (!ind.increaseIndex())
			break;

		for (int idim = 0; idim < dim; idim++)
			x[idim] = x0[idim] + (double) ind.index(idim) * dx[idim];

		if (ind.index(nldim) == 0)
			*stream << std::endl;

	}

}

template<class Grid, typename ft, int fdim> void saveFieldToFile2 (const Grid & g, const Dune::DynamicVector<ft> & F, 
		const Dune::FieldVector<typename Grid::ctype, Grid::dimension> x0,
		const Dune::FieldVector<typename Grid::ctype, Grid::dimension> l,
		const Dune::FieldVector<int, Grid::dimension> N,
		std::string name)
{
	std::ostream * stream = new std::ofstream (name);
	saveField2<Grid,ft,fdim> (g, F, x0, l, N, stream);
	delete stream;
}


#endif
