#ifndef DUNE_INTEGRATE_ENTITY_HH
#define DUNE_INTEGRATE_ENTITY_HH

#include <dune/common/exceptions.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/common/dynvector.hh>

//! compute integral of function over entity with given order
// FIXME: or return FieldVector?
template<class Entity, class Function> 
Dune::FieldVector<typename Function::ftype, Function::fdimension> 
	integrateEntity (const Entity & entity, const Function & func, int p)
{
	// dimension of the entity 
	const int dim = Entity::dimension;
	if (dim != Function::dimension)
		DUNE_THROW(Dune::Exception,"Entity and function's coordinate dimensions disagree");

	// type used for coordinates in the grid
	typedef typename Entity::ctype ctype;

	// get geometry
	const typename Entity::Geometry geometry = entity.geometry();

	// get geometry type
	const Dune::GeometryType gt = geometry.type();

	// get quadrature rule of order p
	const Dune::QuadratureRule<ctype,dim> & rule = Dune::QuadratureRules<ctype,dim>::rule(gt,p);

	// ensure that rule has at least the requested order
	if (rule.order() < p)
		DUNE_THROW(Dune::Exception,"order not available");

	Dune::FieldVector<typename Function::ftype, Function::fdimension> f (0.);

	// compute approximate integral
	for (typename Dune::QuadratureRule<ctype,dim>::const_iterator i=rule.begin(); i!=rule.end(); ++i)
	{
		Dune::FieldVector<typename Function::ftype, Function::fdimension> fval 
//			= func (geometry.global(i->position()), entity);
				= func (i->position(), geometry.global(i->position()), entity);
//				= func (i->position(), entity);


		double weight = i->weight();                 
		double detjac = geometry.integrationElement(i->position());
		for (int i = 0; i < Function::fdimension; i++)
			f[i] += fval[i] * (typename Function::ftype) weight * (typename Function::ftype) detjac;
#ifdef DEBUG
		std::cerr << "weight =" << weight << " detjac =" << detjac << std::endl;
		std::cerr << "fval =" << fval << std::endl;
		std::cerr << "integral = " << f << std::endl;
#endif			
	}

	return f;
}

#endif
