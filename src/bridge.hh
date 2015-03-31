#ifndef CHANNELING_BRIDGE_POTENTIAL_HH
# define CHANNELING_BRIDGE_POTENTIAL_HH

#include"dune/common/fvector.hh" // Field Vector
#include <dune/geometry/quadraturerules.hh>

#ifdef DEBUG
# undef DEBUG
#endif
//#define DEBUG

//#define pow2(x) 	 ((x)*(x))
template <typename ct, int dim, class Function>
class RectBridge
{
	public:

		RectBridge (const Dune::FieldVector<ct, dim> & O, const Dune::FieldVector<ct, dim-1> & w,
			const Function & F, int p)
			: Xo(O), h2(w), F(F), rule(Dune::QuadratureRules<ct,dim-1>::rule(Dune::GeometryType::cube,p))
		{
			// ensure that rule has at least the requested order
			if (rule.order() < p)
				DUNE_THROW(Dune::Exception,"order not available");
			det = 1.;
			for (int i = 0; i < dim-1; i++)
			{
				a[i] = O[i] - w[i];
				k[i] = 2. * w[i];
				det *= k[i];
			}
		};

		 Dune::FieldVector<typename Function::ftype, Function::fdimension> 
		 	operator () (const Dune::FieldVector<ct, dim> & x) const
		{
			Dune::FieldVector<typename Function::ftype, Function::fdimension> f (0.);

			// compute approximate integral
			for (typename Dune::QuadratureRule<ct,dim-1>::const_iterator i=rule.begin(); i!=rule.end(); ++i)
			{
				Dune::FieldVector<ct, dim> y = Xo;
				for (int j = 0; j < dim - 1; j++)
					y[j] = a[j] + k[j] * i->position()[j];
				Dune::FieldVector<typename Function::ftype, Function::fdimension> fval 
					//			= func (geometry.global(i->position()), entity);
					= F(x, y);

				double weight = i->weight();                 
				for (int i = 0; i < Function::fdimension; i++)
					f[i] += fval[i] * (typename Function::ftype) weight * det;
#ifdef DEBUG
				std::cerr << "weight =" << weight << " detjac =" << det << std::endl;
				std::cerr << "fval =" << fval << std::endl;
				std::cerr << "integral = " << f << std::endl;
#endif			
			}

			return f;
		}

	private:

		const Dune::FieldVector<ct, dim> & Xo;    // location of a rectangular bridge 
		const Dune::FieldVector<ct, dim-1> & h2;  // half length from the bridge center 

		// transformation 
		Dune::FieldVector<ct, dim-1> a;   
		Dune::FieldVector<ct, dim-1> k;   
		ct det;

		const Function & F;
		const Dune::QuadratureRule<ct,dim-1> & rule;
};


template <typename ct, int dim, typename ft, int fdim>
class MorseBridge
{
	public:

		MorseBridge (const Dune::FieldVector<ct, dim> & O, const Dune::FieldVector<ct, dim-1> & w,
			ft De, ct a)
			: Xo(O), h2(w), De(De), a(a)
		{
		};

		 Dune::FieldVector<ft, dim*fdim> 
		 	operator () (const Dune::FieldVector<ct, dim> & x) const
		{
			Dune::FieldVector<ft, dim*fdim> f (0.);
			
			for (int i = 0; i < dim - 1; i++)
				if ( (x[i] < Xo[i] - h2[i] - a) || (x[i] > Xo[i] + h2[i] + a) )
					return f;
			// FIXME: no need for fabs()
			double r = fabs(x[dim-1] - Xo[dim-1]);
			// z-component of the force of the first field
			f[dim-1] = 2. * a * De * (exp (-2. * a * r) -  exp (- a * r));
			
			// check if we want to leave the bridge
			// 
#define MOLSIZE			2.  // this is in nm
			if (r < MOLSIZE)
			{
				for (int i = 0; i < dim - 1; i++)
				{	
					// left border
					ct dr = (Xo[i] - h2[i]) - x[i]; 
					if (dr > 0)
					{				
						f[i] = 2. * a * De * (exp (-2. * a * dr) -  exp (- a * dr));
						//std::cerr << "x = " << x << std::endl;
						//std::cerr << "     f (z=" << r << ", dr" << i << "=" << dr << ") = " << f << std::endl;
					}
					// right border
					dr = x[i] - (Xo[i] + h2[i]);
					if (dr > 0)
					{
						f[i] = - 2. * a * De * (exp (-2. * a * dr) -  exp (- a * dr));
						//std::cerr << "x = " << x << std::endl;
						//std::cerr << "     f (z=" << r << ", dr" << i << "=" << dr << ") = " << f << std::endl;
					}

				}
			}
#ifdef DEBUG
			std::cerr << "x = " << x << std::endl;
			std::cerr << "     f ( r=" << r << "; a=" << a << ") = " << f << std::endl;
#endif			
			return f;
		}

	private:

		const Dune::FieldVector<ct, dim> & Xo;    // location of a rectangular bridge 
		const Dune::FieldVector<ct, dim-1> & h2;  // half length from the bridge center 

		ft De;
		ct a;
		};

#endif 

