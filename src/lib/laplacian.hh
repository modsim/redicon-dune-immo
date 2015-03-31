#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

//#define DEBUG
#include "shapeFunctions.hh"
#include "entityFunction.hh"

#ifdef DEBUG
# undef DEBUG
#endif

template<class GridView, typename ft, int fdim>
class entityFunctionLaplacian : public entityP1Function <GridView, ft, fdim>
{
	public:
		entityFunctionLaplacian (int i) {}

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

				Dune::FieldMatrix<typename GridView::ctype, GridView::dimension, GridView::dimension> jacInvTra =
					e.geometry().jacobianInverseTransposed(x);

				for (int i = 0; i < vertexsize; i++)
				{
					Dune::FieldVector<ft,GridView::dimension> grad1;
					jacInvTra.mv(this->basis[i].evaluateGradient(x),grad1);

					for (int j = 0; j < vertexsize; j++)
					{
						Dune::FieldVector<ft,GridView::dimension> grad2;
						jacInvTra.mv(this->basis[j].evaluateGradient(x),grad2);
						for (int a = 0; a < fdim; a++)
							for (int b = 0; b < fdim; b++)
							{
								int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
								if (a == b)
								{
									m[ind] = grad1 * grad2;
									m[ind] *= 1.;
								}
								else 
									m[ind] = 0.;
							}
					}
				}
#ifdef DEBUG
				std::cerr << "m = " << y << std::endl;
#endif			
				return  m;
			}
};

