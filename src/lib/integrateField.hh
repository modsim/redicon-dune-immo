#ifndef DUNE_INTEGRATE_FIELD_HH
# define DUNE_INTEGRATE_FIELD_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

//#define DEBUG
//#define DO_CHECKS
#include "integrate_n.hh"
#include "shapeFunctions.hh"
//#define DEBUG

template<class GridView, typename ft, int fdim>
class entityFunctionIntegrand : public entityP1Function<GridView, ft, fdim>
{
	public:

		enum {fdimension = fdim};

		entityFunctionIntegrand (const GridView & gv, const Dune::DynamicVector<ft> & f) : F(f),  gv(gv), set(gv.indexSet())
		{
				//set = gv.indexSet();
		}

		template<class Entity> Dune::FieldVector<ft, fdim> 
			operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x, 
					const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& X, // global
					const Entity & e) const
			{
#ifdef DEBUG
				std::cerr << "F = " << F << std::endl;
#endif			

				Dune::FieldVector<ft, fdim> V (0.);
				int vertexsize = this->vertexSize(e);
				for (int a = 0; a < fdim; a++)
					for (int i = 0; i < vertexsize; i++)
					{
						int ind2 = set.subIndex(e, i, GridView::dimension)  * fdim + a;
#ifdef DEBUG						
						std::cerr << "a=" << a << " i=" << i << "/" << vertexsize 
							<< ", ind=" << ind2 << " F=" << F[ind2]
							<< ", basis(" << i << ")=" << this->basis[i].evaluateFunction(x) << std::endl;
#endif
						V[a] +=  this->basis[i].evaluateFunction(x) * F[ind2];
					}
#ifdef DEBUG
				std::cerr << "v (@ " << x << ") = " << V << std::endl;
#endif			
				return  V;
			}

	protected:
		const Dune::DynamicVector<ft> & F;
		const GridView & gv;
		const typename GridView::IndexSet & set;
};


template<class Grid, typename ft, int fdim>
	Dune::FieldVector<ft, fdim> 
integrateField (Grid& grid, const Dune::DynamicVector<ft> & F, int p)
{

	// get GridView on leaf grid - type
	typedef typename Grid :: LeafGridView GridView;

	// get GridView instance
	GridView gridView = grid.leafView();

	typedef entityFunctionIntegrand<GridView, ft, fdim> Integrand;
	Integrand f (gridView, F);

	return integrate (grid, f, p);
}
#endif

