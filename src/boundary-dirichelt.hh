#ifndef CHANNELING_IMMO_BOUNDARY_DIRICHLET_HH
# define CHANNELING_IMMO_BOUNDARY_DIRICHLET_HH

#include"dune/grid/common/gridinfo.hh" 
#include"dune/common/fvector.hh" // Field Vector
#include"dune/common/dynmatrix.hh" // Field Vector

//#define DEBUG
#include "setRobin.hh"

#ifdef DEBUG
# undef DEBUG
#endif

/* 
 * Robin BC: d f|n + eta f = u
*/

#define EPS 1.e-8
//#define SOURCE_RADIUS    5.0
//#define SOURCE_CENTER    10.0 // y,z
#define SOURCE_PLANE     2    // on the plane z = 0
#define pow2(x) 	 ((x)*(x))
//#define DEBUG 

template<typename ct, int dim, typename ft, int fdim>
class DirichletBC
{
	public:
		DirichletBC (const Dune::FieldVector<ft, fdim> & u) : k0(u) {};

		DirichletBC (void) : k0 (Dune::FieldVector<double, fdim> (0) ) {};

		//bool hasPoint (const Dune::FieldVector<ct, dim>& x) const {return true;};
		template<class Entity> bool hasEntity (const Entity & e) const 
		{
#ifdef DEBUG		
			std::cerr << "dim= " << Entity::mydimension << std::endl;
			std::cerr << e.geometry().type() << std::endl;
			std::cerr << "center: " << e.geometry().center() << std::endl;
			std::cerr << "corners: " << e.geometry().corners() << std::endl;
			std::cerr << "volume: " << e.geometry().volume() << std::endl;
#endif
			if ( fabs((e.geometry().center())[SOURCE_PLANE]) > EPS ) 
			{
				return true;
			}
			return false;
		}

		Dune::FieldVector<ft, fdim> u (const Dune::FieldVector<ct,dim>& x) const 
		{ 
			
			return k0;
		} ;
	private:
		const Dune::FieldVector<ft, fdim> & k0;

};
#endif // CHANNELING_IMMO_BOUNDARY_HH

