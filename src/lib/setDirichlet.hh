#ifndef DUNE_SET_DIRICHLET_HH
# define DUNE_SET_DIRICHLET_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include"dune/common/dynmatrix.hh"

//#define DEBUG
#include "integrateEntityDynamic.hh"
#include "pointSource.hh"

//#define DEBUG
// WARNING: bc is not implemented, this function sets zero at the Dirichlet Boundary
template<class GridView, typename ftype, int fdim, class DirichletBC> 
bool setDirichlet (const GridView & gv, Dune::DynamicMatrix<ftype> & A, Dune::DynamicVector<ftype> & b, int ifdim, const DirichletBC & bc)
{
	int dim = GridView::dimension;
	typedef typename GridView::ctype ctype;

	// Matrix indeces 
	typedef typename GridView::IndexSet iSet;
	const iSet & set = gv.indexSet();

	if (ifdim >= fdim)
	{
		std::cerr << "setDirichlet(): error: field dimension out of bound (" << ifdim << " > " << fdim << ")"<< std::endl;
		return false;
	}

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
			// determine geometry type of the current element and get the matching reference element
			Dune::GeometryType gt = it->type();
			const Dune::template GenericReferenceElement<ctype,GridView::dimension> &ref =
				Dune::GenericReferenceElements<ctype,GridView::dimension>::general(gt);

				// traverse all vertices the intersection consists of
				for (int i=0; i < ref.size(is->indexInInside(),1,dim); i++)
				{
					// and replace the associated line of A and b with a trivial one
					int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
					indexi = indexi * fdim + ifdim;
//					std::cerr << "index=" << indexi << std::endl;

//					A[indexi] = 0.0; /*@\label{fem:trivialline}@*/
//					A[indexi][indexi] = 1.0;
//					b[indexi] = 0.0;
					
					A[indexi] = 0.0;
					A[indexi][indexi] = 1.0;
					b[indexi] = 0.0;

				}
			}
		}
	}

	return true;
}

#endif

