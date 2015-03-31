#ifndef DUNE_ASSEMBLE_MATRIX_HH
# define DUNE_ASSEMBLE_MATRIX_HH

#include"dune/common/dynmatrix.hh"

//#define DEBUG
#include "integrateEntityDynamic.hh"
#include "pointSource.hh"
#include "updateMatrix.hh"

//#define DEBUG

template<class GridView, class Function, class Source> 
Dune::DynamicMatrix<typename Function::ftype>  assembleMatrix (const GridView & gv, Function & func, int p, std::vector<Source*> slist)
{
	Dune::DynamicMatrix<typename Function::ftype> M = assembleMatrix (gv, func, p);
	updateMatrix<GridView, typename Function::ftype, Function::fdimension, Source> (M, gv, slist);
	return M;
}

template<class GridView, class Function> 
Dune::DynamicMatrix<typename Function::ftype>  assembleMatrix (const GridView & gv, Function & func, int p)
{
	int dim = GridView::dimension;
	int fdim = Function::fdimension;

	Dune::DynamicMatrix<typename Function::ftype> M;
	int N = gv.size (GridView::dimension) * fdim;
	M.resize (N,N);
	M = 0;

	// MAtrix indeces 
	typedef typename GridView::IndexSet iSet;
	const iSet & set = gv.indexSet();

	// get iterator type
	typedef typename GridView :: template Codim<0> :: Iterator LeafIterator;
	for (LeafIterator it = gv.template begin <0>(); it != gv.template end<0>(); ++it)
	{
		// FIXME: Same code in Function::vertexSize(): make final in Cx011
		Dune::GeometryType gt = it->type();
		const Dune::template GenericReferenceElement<typename GridView::ctype, GridView::dimension> &ref =
			Dune::GenericReferenceElements<typename GridView::ctype, GridView::dimension>::general(gt);
		int vertexsize = ref.size(dim);

		Dune::DynamicVector<typename Function::ftype> m =  integrateEntityDynamic (*it, func, p);
#ifdef DO_CHECKS
		if (vertexsize != func.vertexSize(e))
			DUNE_THROW (Dune::Exception, "User provided function gives wrong vertex size");

		int size = vertexsize * fdim;
		size *= size;
		if (m.size() != size)
			DUNE_THROW (Dune::Exception, "User provided function returns DynamicVector of wrong size");
#endif
		for (int i = 0; i < vertexsize; i++)
			for (int j = 0; j < vertexsize; j++)
				for (int a = 0; a < fdim; a++)
					for (int b = 0; b < fdim; b++)
					{
						int ind1 = set.subIndex(*it, i, dim) * fdim  + a;
						int ind2 = set.subIndex(*it, j, dim) * fdim  + b;
						int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;

						M[ind1][ind2] += m[ind];
					}
	}

	//updateMatrix<GridView, typename Function::ftype, Function::fdimension, Source> (M, gv, slist);

/*
//#define DO_CHECKS
	for (typename std::vector<Source*>::iterator ps = slist.begin(); ps != slist.end(); ++ps)
	{	
		Source * is = *ps;
		Dune::DynamicVector<typename Function::ftype> m =  (*is) ();
		int vertexsize = is->entityVertexSize ();

		unsigned int size = vertexsize * fdim;
		size *= size;
#ifdef DO_CHECKS
		if (m.size() != size)
			DUNE_THROW (Dune::Exception, "Point Sources: User provided function returns DynamicVector of wrong size");
#endif
		for (int i = 0; i < vertexsize; i++)
			for (int j = 0; j < vertexsize; j++)
				for (int a = 0; a < fdim; a++)
					for (int b = 0; b < fdim; b++)
					{
						int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;

						int ind1 = is->entityIndex (i) * fdim  + a;
						int ind2 = is->entityIndex(j) * fdim  + b;

						M[ind1][ind2] += m[ind];
						std::cerr << "M[" << ind1 << "][" << ind2 
							<< "]+= " << m[ind] << std::endl;

					}
	}
*/
	return M;
}
#undef DO_CHECKS
#endif

