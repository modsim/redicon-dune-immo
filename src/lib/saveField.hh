#ifndef  DUNE_SAVE_FIELD_HH
# define DUNE_SAVE_FIELD_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <iostream>
#include <fstream>

#include "dune/grid/common/mcmgmapper.hh" 
#include "dune/common/dynvector.hh" 

template<int dim> struct VertexLayout 
{
	bool contains (Dune::GeometryType gt)
	{
		if (gt.dim() == 0)
			return true;

		return false;
	}
};

template<class Grid, typename ft, int fdim> void saveField (const Grid & g, const Dune::DynamicVector<ft> & f, std::ostream * stream)
{
	const int dim = Grid::dimension;
	typedef typename Grid::LeafGridView GridView;
	const GridView & gv = g.leafView();
	typedef typename GridView::template Codim<dim>::Iterator Iterator; 
	Dune::LeafMultipleCodimMultipleGeomTypeMapper<Grid, VertexLayout> mapper (g);

	//std::cerr << "f=" << f << std::endl;

	for (Iterator it = gv.template begin<dim>(); it!=gv.template end<dim>(); ++it)
	{
		*stream << (*it).geometry().corner(0);

		int index = mapper.map (*it);
		//std::cerr << "index=" << index ;

		for (int a = 0; a < fdim; ++a)
		{
			int ind = index * fdim + a;
	//		std::cerr << "F[" << ind << "]=" << f[ind] << std::endl;
			*stream << "	" << f[ind];

		}
		*stream << std::endl;
	}
}

template<class Grid, typename ft, int fdim> void saveFieldToFile (const Grid & g, const Dune::DynamicVector<ft> & f, std::string name)
{
	//std::cerr << "name=" << name << std::endl;
	std::ostream * stream = new std::ofstream (name);
	saveField<Grid,ft,fdim> (g, f, stream);
	delete stream;
}


#endif
