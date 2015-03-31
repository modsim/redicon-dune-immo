#ifndef DUNE_SURFACE_AREA_HH
# define DUNE_SURFACE_AREA_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include"dune/common/dynmatrix.hh"

//#define DEBUG
#include "integrateEntityDynamic.hh"
#include "pointSource.hh"

//#define DEBUG
// WARNING: bc is not implemented, this function sets zero at the Dirichlet Boundary
template<class GridView, class DirichletBC> 
typename GridView::ctype calcSurfaceArea (const GridView & gv, const DirichletBC & bc)
{
	typedef typename GridView::ctype ctype;
	ctype area = 0.;;

	// Matrix indeces 
	//typedef typename GridView::IndexSet iSet;
	//const iSet & set = gv.indexSet();

	// get iterator type
	typedef typename GridView :: template Codim<0> :: Iterator LeafIterator;
	for (LeafIterator it = gv.template begin <0>(); it != gv.template end<0>(); ++it)
	{
		typedef typename GridView::IntersectionIterator IntersectionIterator;

		const IntersectionIterator isend = gv.iend(*it);
		for (IntersectionIterator is = gv.ibegin(*it) ; is != isend ; ++is)
		{
			// check whether current intersection is on the correct boundary
			if ( is->boundary()  && bc.hasEntity (*is) )
			{
				area += (*is).geometry().volume();
			}
		}
	}

	return area;
}

#endif

