#ifndef DUNE_SOURCE_FUNCTIONS_HH
# define DUNE_SOURCE_FUNCTIONS_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif
#include <iostream>
#include"dune/common/mpihelper.hh" // An initializer of MPI
#include"dune/common/exceptions.hh" // We use exceptions
#include"dune/grid/onedgrid.hh" //  1D Grid
#include"dune/grid/common/gridinfo.hh" //  1D Grid
#include"dune/common/fvector.hh" // Field Vector
#include"dune/common/dynmatrix.hh" // Field Vector

//#define DEBUG
#include "shapeFunctions.hh"

#ifdef DEBUG
# undef DEBUG
#endif
//#define DEBUG

template <class GridView,  typename ft, int fdim> 
class SourceFunctionMatrix : public entityP1Function <GridView, ft, fdim>
{
	public: 

		enum {fdimension = fdim};
		typedef ft ftype;

		SourceFunctionMatrix (const Dune::FieldMatrix<ft, fdim, fdim> & Q) : q(Q) {};
		
		template<class Entity> int size (const Entity & e) const
		{
				int size = fdim * this->vertexSize (e);
				return size * size;
		}

		Dune::FieldMatrix<ft, fdim, fdim> Q (int key) const {return q;};

		template<class Entity> Dune::DynamicVector<ft> 
		operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x, const Entity & e, int key) const
		{
			Dune::DynamicVector<ft> M;
			int vertexsize = this->vertexSize(e);
			M.resize (size(e));
			M = 0;

			// get geometry
			// const typename Entity::Geometry geometry = e.geometry();
			//const Dune::FieldVector<typename GridView::ctype,GridView::dimension> & y = geometry.local(x);
			const Dune::FieldVector<typename GridView::ctype,GridView::dimension> & y = x;

			for (int i = 0; i < vertexsize; i++)
				for (int j = 0; j < vertexsize; j++)
					for (int a = 0; a < fdim; a++)
						for (int b = 0; b < fdim; b++)
						{
							int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
							M[ind] = this->q[a][b] * this->basis[i].evaluateFunction(y) * this->basis[j].evaluateFunction(y);

#ifdef DEBUG						
							std::cerr << "x=" << x 
								<< " basis(" << i << ")=" << this->basis[i].evaluateFunction(y)
								<< " basis(" << j << ")=" << this->basis[j].evaluateFunction(y) 
								<< std::endl;
#endif

						}
#ifdef DEBUG
			std::cerr << "m = " << M << std::endl;
#endif		
			return  M;
		}

	private:
		const Dune::FieldMatrix<ft, fdim, fdim> & q;
};

template <class GridView,  typename ft, int fdim> 
class SourceFunctionVector : public entityP1Function <GridView, ft, fdim>
{
	public: 
		SourceFunctionVector (const Dune::FieldVector<ft, fdim> & Q) : q(Q) {};

		
		template<class Entity> int size (const Entity & e) const
		{
				return fdim * this->vertexSize (e);
		}

		template<class Entity> Dune::DynamicVector<ft> 
		operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x, const Entity & e, int key) const
		{
			Dune::DynamicVector<ft> V ;
			int vertexsize = this->vertexSize(e);
			V.resize (size(e));
			V = 0;

			// get geometry
			const typename Entity::Geometry geometry = e.geometry();
			//const Dune::FieldVector<typename GridView::ctype,GridView::dimension> & y = geometry.local(x);
			const Dune::FieldVector<typename GridView::ctype,GridView::dimension> & y = x;

			for (int i = 0; i < vertexsize; i++)
				for (int a = 0; a < fdim; a++)
				{
					int ind2 = i * fdim + a;
					V[ind2] = this->q[a] * this->basis[i].evaluateFunction(y);
#ifdef DEBUG
					std::cerr << "i=" << i << " a=" << a
					<< " ind=" << ind2
					<< " basis=" << this->basis[i].evaluateFunction(y)
					<< " v= " << V[ind2] << std::endl;
#endif
				}
#ifdef DEBUG
			std::cerr << "sfunc(): v = " << V << std::endl;
#endif		
			return  V;
		}
	private:
		const Dune::FieldVector<ft, fdim> & q;
};

#endif
