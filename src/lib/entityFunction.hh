#ifndef DUNE_ENTITY_FUNCTION_HH
# define DUNE_ENTITY_FUNCTION_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include "shapeFunctions.hh"
#include"dune/common/dynmatrix.hh" // Field Vector

template<class GridView, typename ft, int fdim>
class entityP1Function
{
	public:
		entityP1Function () 
		{
			basis = 
			P1ShapeFunctionSet<typename GridView::ctype, typename GridView::ctype, GridView::dimension>::instance();
		}

		template<class Entity> const int vertexSize (const Entity & e) const
		{ 
			Dune::GeometryType gt = e.type();
			const Dune::template GenericReferenceElement<typename GridView::ctype, GridView::dimension> &ref =
				Dune::GenericReferenceElements<typename GridView::ctype, GridView::dimension>::general(gt);
			return ref.size(GridView::dimension);
		};

		P1ShapeFunctionSet<typename GridView::ctype, typename GridView::ctype, GridView::dimension> basis;

		enum {dimension = GridView::dimension};
		enum {fdimension = fdim};
		typedef typename GridView::ctype ctype;
		typedef ft ftype;


};
#endif

