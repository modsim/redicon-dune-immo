#ifndef DUNE_SET_ROBIN_HH
# define DUNE_SET_ROBIN_HH

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include"dune/common/dynmatrix.hh"

//#define DEBUG
#include "integrateEntityDynamic.hh"
#include "entityFunction2.hh"

//#define DEBUG

/* 
 * Robin BC: d f|n + eta f = u
*/

//#define DEBUG

template<typename ct, int dim, typename ft, int fdim, class RobinBC>
class entityFuncRobinMatrix : public entityP1Function2 <ct, dim, ft, fdim>
{
	public:
		entityFuncRobinMatrix (const RobinBC & BC) : BC(BC) {};

		template<class Entity> int size (const Entity & e) const { const int a = fdim * this->vertexSize (e); return a*a;} 
		template<class Entity> Dune::DynamicVector<ft> 
			operator() (const Dune::FieldVector<ct, dim> & x, 
					const Dune::FieldVector<ct, dim+1> & X,
					const Entity & e) const
			{
				Dune::DynamicVector<ft> m ;
				int vertexsize = this->vertexSize (e) ;
				m.resize (size(e));
				m = 0;
#ifdef DEBUG
				std::cerr << "entityFuncRobinMatrix(): x = " << x << std::endl;
#endif			
				const typename Entity::Geometry geometry = e.geometry();
				double detjac = geometry.integrationElement (x);
				//double detjac = 1.;
#ifdef DEBUG
				std::cerr << "entityFuncRobinMatrix(): detjac= " << detjac << std::endl;
#endif			

				for (int i = 0; i < vertexsize; i++)
				{
					//Dune::FieldVector<ft,fdim> N1 = this->basis[i].evaluateFunction(x);
					Dune::FieldVector<double,1> N1 = this->basis[i].evaluateFunction(x);

					for (int j = 0; j < vertexsize; j++)
					{
						//Dune::FieldVector<ft,fdim> N2 = this->basis[j].evaluateFunction(x) ;
						Dune::FieldVector<double,1> N2 = this->basis[j].evaluateFunction(x) ;
						for (int a = 0; a < fdim; a++)
							for (int b = 0; b < fdim; b++)
							{
								int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
								 m[ind] =  BC.eta(X)[a][b] * N1 * N2 / detjac / detjac;
								//if (a == b)
								//{
								 //	m[ind] =  (BC.eta(x))[a] * N1 * N2 / detjac;
								//}
								//else 
								//	m[ind] = 0.;

							}
					}
				}
#ifdef DEBUG
				std::cerr << "G = " << m << std::endl;
#endif			
				return  m;
			}
	private:
		const RobinBC & BC;
};


template<typename ct, int dim, typename ft, int fdim, class RobinBC>
class entityFuncRobinVector : public entityP1Function2 <ct, dim, ft, fdim>
{
	public:
		entityFuncRobinVector (const RobinBC & BC) : BC(BC) {}

		template<class Entity> int size (const Entity & e) const { const int a = fdim * this->vertexSize (e); return a;} 

		template<class Entity> Dune::DynamicVector<ft> 
			operator() (const Dune::FieldVector<ct, dim>& x, 
				const Dune::FieldVector<ct, dim+1>& X,
				const Entity & e) const
			{
				Dune::DynamicVector<ft> v ;
				int vertexsize = this->vertexSize (e) ;
				v.resize (size(e));
				v = 0;
#ifdef DEBUG
				std::cerr << "entityFuncRobinVector(): x = " << x << std::endl;
#endif			
				const typename Entity::Geometry geometry = e.geometry();
				//double detjac = geometry.integrationElement (x);
				//double detjac = 1.;

				for (int i = 0; i < vertexsize; i++)
				{
					//Dune::FieldVector<ft, fdim> N1 = this->basis[i].evaluateFunction(x);
					Dune::FieldVector<double, 1> N1 = this->basis[i].evaluateFunction(x);
					for (int a = 0; a < fdim; a++)
					{
						int ind2 = i * fdim + a;
						//v[ind2] =  BC.u(x) * N1  / detjac;
//						std::cerr << "u(" << a << ") = " << (BC.u(x))[a] << std::endl;
						v[ind2] =  (BC.u(X))[a] * N1;
					}
				}
#ifdef DEBUG
				std::cerr << "V = " << v << std::endl;
#endif			
				return  v;
			}
	private:
		const RobinBC & BC;
};

// FIXME: make RobinBC a template parameter
template<class GridView, typename ftype, int fdim, class RobinBC> 
bool setRobin (const GridView & gv, Dune::DynamicMatrix<ftype> & A, Dune::DynamicVector<ftype> & B, 
//		const RobinBC<typename GridView::ctype, GridView::dimension - 1, ftype, fdim> & bc)
		const RobinBC & bc)
{
	// FIXME: make an argument
	int p = 6;

	const int dim = GridView::dimension ;
	typedef typename GridView::ctype ctype;

	entityFuncRobinMatrix<ctype, dim - 1, double, fdim, RobinBC> Rm (bc);
	entityFuncRobinVector<ctype, dim - 1 , double, fdim, RobinBC> Rv (bc);

	// Matrix indeces 
	typedef typename GridView::IndexSet iSet;
	const iSet & set = gv.indexSet();

#ifdef DEBUG
	std::cerr << "Starting loop " << std::endl;
#endif
	// get iterator type
	typedef typename GridView :: template Codim<0> :: Iterator LeafIterator;
	for (LeafIterator it = gv.template begin <0>(); it != gv.template end<0>(); ++it)
	{
		typedef typename GridView::IntersectionIterator IntersectionIterator;
		for (IntersectionIterator is = gv.ibegin(*it) ; is != gv.iend(*it) ; ++is)
		{
			// determine geometry type of the current element and get the matching reference element

			// check whether current intersection is on the boundary
			// and has Robin BC
			if ( is->boundary() && bc.hasEntity (*is))
			{
				Dune::DynamicVector<ftype> m =  integrateEntityDynamic (*is, Rm, p);
				Dune::DynamicVector<ftype> v =  integrateEntityDynamic (*is, Rv, p);
#ifdef DEBUG
				std::cerr << "*** Boundary:" << std::endl;
				std::cerr << "m=" << m << std::endl;
				std::cerr << "v=" << v << std::endl;
#endif	

				Dune::GeometryType gt = it->type();
				const Dune::template GenericReferenceElement<ctype, dim> &ref =
					Dune::GenericReferenceElements<ctype, dim>::general(gt);

				// traverse all vertices the intersection consists of
				//  FIXME: matrix is symmetric -- optimize
				int vertexsize = ref.size(is->indexInInside(),1,dim);
				for (int i=0; i < vertexsize; i++)
				{
					int indexi = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,i,dim),dim);
#ifdef DEBUG					
					std::cerr << "boundary indexi=" << indexi;
#endif
					for (int a = 0; a < fdim; a++)
					{
						int indv = i * fdim + a;
						B[indexi * fdim  + a] += v[indv];
#ifdef DEBUG						
						std::cerr << " field index a=" << a << " adding v=" <<  v[indv] << std::endl;
#endif					
						for (int j = 0; j < vertexsize; j++)
						{
							int indexj = set.subIndex(*it,ref.subEntity(is->indexInInside(),1,j,dim),dim);
#ifdef DEBUG							
							std::cerr << "boundary indexj=" << indexj;
#endif
							for (int b = 0; b < fdim; b++)
							{
								int indm = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;
								A[indexi * fdim  + a][indexj * fdim  + b] += m[indm];
#ifdef DEBUG								
								std::cerr << "field index b=" << b << " m=" <<  m[indm] << std::endl;
#endif								
							}
						}
					}
				}
			}

		}
	}
#ifdef DEBUG	
	std::cerr << "B = " << B << std::endl;
#endif

	return true;
}

#ifdef DEBUG
#  undef DEBUG
#endif

#endif

