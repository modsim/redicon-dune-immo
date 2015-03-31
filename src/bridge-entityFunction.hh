#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

//#define DEBUG
#include "shapeFunctions.hh"
#include "entityFunction.hh"

#ifdef DEBUG
# undef DEBUG
#endif

//
// This is the entity function for the current
//    Ja = Da (r) \nabla Phi_a + Phi_a x Fa(r)
// where Da is the diffusion coefficient, and Fa force

template<class GridView, typename ft, int fdim, class Diffusion, class Force>
class entityFunctionForce : public entityP1Function <GridView, ft, fdim>
{
	public:
		entityFunctionForce (const Diffusion & D, const Force & F) : D(D), F(F) {}

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

				const Dune::FieldMatrix<ft,GridView::dimension*fdim, GridView::dimension*fdim> & diff = D(X);
				const Dune::FieldVector<ft,GridView::dimension*fdim> & force = F(X);

				Dune::FieldMatrix<typename GridView::ctype, GridView::dimension, GridView::dimension> jacInvTra =
					e.geometry().jacobianInverseTransposed(x);

				for (int i = 0; i < vertexsize; i++)
				{
					Dune::FieldVector<ft,GridView::dimension> grad1;
					jacInvTra.mv(this->basis[i].evaluateGradient(x),grad1);
					//Dune::FieldVector<double, 1> N1 = this->basis[i].evaluateFunction(x);

					for (int j = 0; j < vertexsize; j++)
					{
						Dune::FieldVector<ft,GridView::dimension> grad2;
						jacInvTra.mv(this->basis[j].evaluateGradient(x),grad2);
						Dune::FieldVector<double, 1> N2 = this->basis[j].evaluateFunction(x);

						// totaly general (but computationaly not too efficient)
						for (int a = 0; a < fdim; a++)
							for (int b = 0; b < fdim; b++)
							{
								int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
								m[ind] = 0;
								for (int k = 0; k < GridView::dimension; k++)
								{
									for (int n = 0; n < GridView::dimension; n++)
									{
										m[ind] += diff[GridView::dimension*a + k][GridView::dimension*b + n] 
											* grad1[k] * grad2[n];
									}
									if (a == b )
										m[ind] += force[GridView::dimension*a + k] * grad1[k] * N2;
								}
							}
					}
				}
#ifdef DEBUG
				std::cerr << "m = " << y << std::endl;
#endif			
				return  m;
			}
		private:
			const Diffusion & D;
			const Force & F;
};

