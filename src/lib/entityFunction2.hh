#ifndef DUNE_ENTITY_FUNCTION_2_HH
# define DUNE_ENTITY_FUNCTION_2_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include "shapeFunctions.hh"
#include"dune/common/dynmatrix.hh" // Field Vector

template<typename ct, int dim, typename ft, int fdim>
class entityP1Function2
{
	public:
		entityP1Function2 () 
		{
			basis = 
			P1ShapeFunctionSet<ct, ct, dim>::instance();
		}

		template<class Entity> const int vertexSize (const Entity & e) const
		{ 
			Dune::GeometryType gt = e.type();
			const Dune::template GenericReferenceElement<ct, dim> &ref =
				Dune::GenericReferenceElements<ct, dim>::general(gt);
			return ref.size(dim);
		};

		P1ShapeFunctionSet<ct, ct, dim> basis;

		enum {dimension = dim};
		enum {fdimension = fdim};
		typedef ct ctype;
		typedef ft ftype;


};
#endif

