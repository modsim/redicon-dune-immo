#ifndef DUNE_EVAL_FIELD_HH
# define DUNE_EVAL_FIELD_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include"dune/common/fvector.hh" // Field Vector
#include "entityFunction.hh"

//#define DEBUG
#include "insideEntity.hh"

//#define DEBUG
template<class GridView, typename ft, int fdim>
class evalField : public entityP1Function<GridView, ft, fdim>
{
	public:
		evalField (const GridView & gv, const Dune::DynamicVector<ft> & f) : F(f),  gv(gv), set(gv.indexSet()) {};
	
		//template<class Entity> 
		Dune::FieldVector<ft, fdim> 
			operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x) const
			{
				typedef typename GridView ::template Codim<0>::Entity Entity;
			//	typedef typename GridView :: template Codim<0> :: Iterator Entity;
				Entity * e = NULL;
				Dune::FieldVector<typename GridView::ctype,GridView::dimension> local;
				for (typename GridView :: template Codim<0> :: Iterator it = gv. template begin <0>(); it != gv. template end<0>(); ++it)
				{
#ifdef DEBUG				
					std::cerr << "in the loop: x=" << x << std::endl;
#endif					
					if (insideEntity<typename GridView :: template Codim<0> :: Iterator, typename GridView::ctype, GridView::dimension> (x, it))
					{
						e = &(*it);
						local = it->geometry().local(x);
#ifdef DEBUG						
						std::cerr << "found!" << std::endl;
#endif

						break;
					}
				}

				Dune::FieldVector<ft, fdim> V (0.);

				if (!e)
				{
//					DUNE_THROW (Dune::Exception, "Cannot determine my entity in grid, probably not in grid");
//					DUNE_THROW (Dune::Exception, "Cannot determine my entity in grid, probably not in grid");
					std::cerr << "evalField(): point " << x << " not in mesh?" << std::endl;
					return V;
				}
#ifdef DEBUG
				std::cerr << "F = " << F << std::endl;
#endif			

				int vertexsize = this->vertexSize(*e);
				for (int i = 0; i < vertexsize; i++)
				{
					for (int a = 0; a < fdim; a++)
					{
						int ind2 = set.subIndex(*e, i, GridView::dimension)  * fdim + a;
#ifdef DEBUG						
						std::cerr << "a=" << a << " i=" << i << "/" << vertexsize 
							<< ", ind=" << ind2 << " F=" << F[ind2]
							<< ", basis(" << i << ")=" << this->basis[i].evaluateFunction(local) << std::endl;
#endif
						V[a] +=  this->basis[i].evaluateFunction(local) * F[ind2];
					}
#ifdef DEBUG
					std::cerr << "v (@ " << x << ") = " << V << std::endl;
#endif			
				}
				return  V;

			}

	protected:
		const Dune::DynamicVector<ft> & F;
		const GridView & gv;
		const typename GridView::IndexSet & set;
};

#ifdef DEBUG
# undef DEBUG
#endif

#endif
