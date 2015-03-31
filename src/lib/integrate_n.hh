#ifndef DUNE_INTEGRATE_N_HH
#define DUNE_INTEGRATE_N_HH

// Dune includes
#include"config.h"           // file constructed by ./configure script
#include"integrateEntity.hh"

//! uniform refinement test
template<class Grid, class Function> 
Dune::FieldVector<typename Function::ftype, Function::fdimension> integrate (Grid& grid, Function & func, int p) 
{

	// get GridView on leaf grid - type
	typedef typename Grid :: LeafGridView GridView;

	// get GridView instance
	GridView gridView = grid.leafView();

	// get iterator type
	typedef typename GridView :: template Codim<0> :: Iterator LeafIterator;

	Dune::FieldVector<typename Function::ftype, Function::fdimension> val (0.);

	// compute integral with some order
	for (LeafIterator it = gridView.template begin<0>(); it != gridView.template end<0>(); ++it)
	{
		Dune::FieldVector<typename Function::ftype, Function::fdimension> value = integrateEntity(*it, func, p);
		val += value;
//		for (int i = 0 ; i < Function::fdimension::; i++)
//			val[i] += value[i];
	}

	return val;
}

#endif
