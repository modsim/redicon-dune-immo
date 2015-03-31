#ifndef CHANNELING_IMMO_BOUNDARY_HH
# define CHANNELING_IMMO_BOUNDARY_HH

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

#define EPS_ROBIN 1.e-4
//#define SOURCE_RADIUS    5.0
//#define SOURCE_CENTER    10.0 // y,z
#define SOURCE_PLANE     2    // on the plane z = 0
#define pow2(x) 	 ((x)*(x))
//#define DEBUG 

template<typename ct, int dim, typename ft, int fdim, class Force>
class RobinBCForce
{
	public:

		RobinBCForce (const Dune::FieldVector<ct, dim+1> & Xo, double Ro, 
				const Dune::FieldVector<ft, fdim> & u, const Dune::FieldMatrix<ft, fdim, fdim> & eta,
				const Force & F) 
			: Xo(Xo), Ro(Ro), k0(u), k1(eta), F(F)
		{};

		RobinBCForce (const Dune::FieldVector<ct, dim+1> & Xo, double Ro, 
				const Dune::FieldVector<ft, fdim> & u,
				const Force & F) 
			: Xo(Xo), Ro(Ro), k0(u), k1(Dune::FieldMatrix<double, fdim,fdim> (0)), F(F)
		{};

		RobinBCForce (const Dune::FieldVector<ct, dim+1> & Xo, double Ro, 
				const Dune::FieldMatrix<ft, fdim, fdim> & eta,
				const Force & F) 
			: Xo(Xo), Ro(Ro), k0(Dune::FieldVector<double, fdim> (0)), k1(eta), F(F)
		{};

		bool hasPoint (const Dune::FieldVector<ct, dim>& x) const {return true;};
		template<class Entity> bool hasEntity (const Entity & e) const 
		{
#ifdef DEBUG		
//			std::cerr << "dim= " << Entity::mydimension << std::endl;
//			std::cerr << e.geometry().type() << std::endl;
//			std::cerr << "center: " << e.geometry().center() << std::endl;
//			std::cerr << "corners: " << e.geometry().corners() << std::endl;
#endif

			if ( fabs((e.geometry().center())[SOURCE_PLANE] - Xo[SOURCE_PLANE]) > EPS_ROBIN ) 
			{
#ifdef DEBUG			
//				std::cerr << "Wrong surface, FALSE "  << std::endl;
#endif				
				return false;
			}
#ifdef DEBUG			
			else		
//				std::cerr << "Surface OK"  << std::endl;
#endif

			for (int i = 0; i < e.geometry().corners(); i++)
			{
#ifdef DEBUG				
				std::cerr << "Corner " << i <<  " : " << e.geometry().corner(i) << std::endl;
#endif
				for (int j = 0; j < dim + 1; j++) // coordinates are in surface dimension (dim) + 1
				{
					ct r = e.geometry().corner(i)[j] - Xo[j];
					if ( (j != SOURCE_PLANE) && ( fabs(r) > Ro + EPS_ROBIN) )
					{
#ifdef DEBUG				
						std::cerr << " dim " << j << " : " << e.geometry().corner(i) 
							<< " outside robin circle " << Xo
							<< ": " << r << " > " << Ro
							<< std::endl;	
#endif
						return false;
					}
#ifdef	DEBUG
					else
						std::cerr << " dim " << j << " OK"  << std::endl;
#endif					
				}
			}
#ifdef DEBUG
			std::cerr << "BOUNDARY: " << std::endl;;
			for (int i = 0; i < e.geometry().corners(); i++)
				std::cerr << "	-- corner " << i << ": " << e.geometry().corner(i) << std::endl;
#endif
			return true;

		};

		Dune::FieldVector<ft, fdim> u (const Dune::FieldVector<ct,dim+1>& x) const 
		{ 
			return k0;
		} ;

		Dune::FieldMatrix<ft, fdim, fdim> eta (const Dune::FieldVector<ct, dim+1>& x) const 
		{
			Dune::FieldVector<ft, fdim * (dim+1)> forceVal = F(x);
			//std::cerr << "robinBC(): force=" << forceVal << std::endl;

			Dune::FieldMatrix<ft, fdim, fdim> etaVal = k1;
			for (int a = 0; a < fdim; a++)
				etaVal[a][a] -= forceVal[a * (dim +1) + dim]; // z-component (dim is here 2 in 3D space)

			return etaVal;
		};

	private:
		const Dune::FieldVector<ct, dim+1> & Xo;
		const double Ro;
		const Dune::FieldVector<ft, fdim> & k0;
		const Dune::FieldMatrix<ft, fdim, fdim> & k1;
		const Force & F;

};
#undef EPS_ROBIN
#endif // CHANNELING_IMMO_BOUNDARY_HH

