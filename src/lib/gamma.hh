#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

//#define DEBUG
#include "shapeFunctions.hh"
#include "entityFunction.hh"

#ifdef DEBUG
# undef DEBUG
#endif
//#define DEBUG
template<class GridView, typename ft, int fdim>
class entityFunctionGamma : public entityP1Function <GridView, ft, fdim>
{
	public:
		entityFunctionGamma (int i) {}

		template<class Entity> int size (const Entity & e) const { const int a = fdim * this->vertexSize (e); return a*a;} 

		template<class Entity> Dune::DynamicVector<ft> 
			operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x, // local
					const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& X, // global
					const Entity & e) const
			{
				Dune::DynamicVector<ft> m ;
				int vertexsize = this->vertexSize (e) ;
				m.resize (size(e));
				m = 0;
#ifdef DEBUG
				std::cerr << "x = " << x << std::endl;
#endif			

	//			const typename Entity::Geometry geometry = e.geometry();
	//			double detjac = geometry.integrationElement (x);

				for (int i = 0; i < vertexsize; i++)
				{
					Dune::FieldVector<ft,1> N1 = this->basis[i].evaluateFunction(x);

					for (int j = 0; j < vertexsize; j++)
					{
						Dune::FieldVector<ft,1> N2 = this->basis[j].evaluateFunction(x) ;
						for (int a = 0; a < fdim; a++)
							for (int b = 0; b < fdim; b++)
							{
								int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
								if (a == b)
								//m[ind] =  N1 * N2 / detjac / detjac;
									m[ind] =  N1 * N2 ;
								else 
									m[ind] = 0.;

							}
					}
				}
#ifdef DEBUG
				std::cerr << "G = " << m << std::endl;
#endif			
				return  m;
			}
};
//#undef DEBUG

