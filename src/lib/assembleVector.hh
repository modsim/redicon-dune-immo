#ifndef DUNE_ASSEMBLE_VECTOR_HH
# define DUNE_ASSEMBLE_VECTOR_HH

//#define DEBUG
#define DO_CHECKS
#include "integrateEntityDynamic.hh"

template<class GridView, class Function>
Dune::DynamicVector<typename Function::ftype> assembleVector (const GridView & gv, Function & func, int p)
{
	int dim = GridView::dimension;
	int fdim = Function::fdimension;

	Dune::DynamicVector<typename Function::ftype> V;
	int N = gv.size (dim)*fdim;
	V.resize (N);

	// Matrix indeces 
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

		Dune::DynamicVector<typename Function::ftype> v =  integrateEntityDynamic (*it, func, p);
#ifdef DO_CHECKS

		if (vertexsize != func.vertexSize(*it))
			DUNE_THROW (Dune::Exception, "User provided function gives wrong vertex size");

		unsigned int size = vertexsize * fdim;
		if (v.size() != size)
			DUNE_THROW (Dune::Exception, "User provided function returns DynamicVector of wrong size");
#endif
		for (int i = 0; i < vertexsize; i++)
			for (int a = 0; a < fdim; a++)
				V[set.subIndex(*it, i, dim) * fdim + a]  += v[i * fdim + a];
	}

	return V;
}


// FIXME add this return type and revertVector()
//std::vector<class sourceVecorElement*>

template<class GridView, typename ftype, int fdim, class Source> 
void updateVector (Dune::DynamicVector<ftype> & V, const GridView & gv, std::vector<Source*> slist)
{

#define DO_CHECKS
	for (typename std::vector<Source*>::iterator ps = slist.begin(); ps != slist.end(); ++ps)
	{	
		Source * is = *ps;
		Dune::DynamicVector<ftype> v =  (*is) ();
		int vertexsize = is->entityVertexSize ();

#ifdef DO_CHECKS
		unsigned int size = vertexsize * fdim;
		if (v.size() != size)
			DUNE_THROW (Dune::Exception, "Point Sources: User provided function returns DynamicVector of wrong size");
#endif
		for (int i = 0; i < vertexsize; i++)
			for (int a = 0; a < fdim; a++)
			{
#ifdef DEBUG
				int ind = i * fdim + a;
				std::cerr << "i=" << i << " a=" << a
				<< " ind=" << ind << " (global_ind=" << is->entityIndex (i) * fdim  + a
				<< ") v= " << v[i * fdim + a] << std::endl;
#endif				
				V[is->entityIndex (i) * fdim  + a] += v[i * fdim + a];
			}
	}
	return;
}


template<class GridView, class Function, class Source>
Dune::DynamicVector<typename Function::ftype> assembleVector (const GridView & gv, Function & func, int p, std::vector<Source*> slist)
{
	Dune::DynamicVector<typename Function::ftype> V = assembleVector (gv, func, p);
	updateVector<GridView, Function::ftype, Function::fdimension. Source> (V, gv, slist);
	return V;
}


#endif

