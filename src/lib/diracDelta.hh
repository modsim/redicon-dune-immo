#ifndef DUNE_DIRACT_DELTA_HH
# define DUNE_DIRAC_DELTA_HH

#include"dune/common/fvector.hh" // Field Vector
//#define DEBUG
#include "insideEntity.hh"

template<typename ct, int dim> 
class diracDelta
{
	public:

		diracDelta (const Dune::FieldVector<ct, dim> & x)
//		diracDelta (Dune::FieldVector<ct, dim> x)
		{
			global_ = x;
		}

		enum {dimension = dim};
		typedef ct ctype;

		void resetPosition (const Dune::FieldVector<ct, dim> & x) {global_ = x;};

		const Dune::FieldVector<ct, dim> x () const {return global_;}
		const Dune::FieldVector<ct, dim> global () const {return global_;}

		template<class Entity> const bool inside (const Entity & e) const
		{ 
//			std::cerr << "deta(x-" << global_ << ")" << std::endl;
			return insideEntity<Entity, ctype, dimension> (global_, e);
		};

		private:
			Dune::FieldVector<ct, dim> global_;
};
#endif

