#ifndef DUNE_INTEGRATE_ENTITY_DYNAMIC_HH
#define DUNE_INTEGRATE_ENTITY_DYNAMIC_HH

#include<dune/common/exceptions.hh>
#include<dune/geometry/quadraturerules.hh>
#include"dune/common/dynvector.hh"

//#define DEBUG

//! compute integral of function over entity with given order
// FIXME: or return FieldVector?
template<class Entity, class Function> 
Dune::DynamicVector<typename Function::ftype> 
	integrateEntityDynamic (const Entity & entity, const Function & func, int p)
{
	// dimension of the entity 
//	const int dim = Entity::dimension;
	const int dim = Entity::mydimension;
#ifdef DEBUG
	std::cerr << "integrateDynamicEntity(): dim= " << Entity::dimension << std::endl;
#endif
	if (dim != Function::dimension)
		DUNE_THROW(Dune::Exception,"Entity and function's coordinate dimensions do not agree");

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

	Dune::DynamicVector<typename Function::ftype> f;
	int size = func.size(entity);
	f.resize (size);

	// compute approximate integral
	for (typename Dune::QuadratureRule<ctype,dim>::const_iterator r=rule.begin(); r!=rule.end(); ++r)
	{
		Dune::DynamicVector<typename Function::ftype> fval
//			= func (geometry.global(i->position()), entity);
			= func (r->position(), geometry.global(r->position()), entity);
//			= func (r->position(), entity);
#ifdef DEBUG			
		std::cerr << "r =" << r->position() 
			<< " -> global = " << geometry.global(r->position()) << std::endl;
#endif
		double weight = r->weight();                 
		double detjac = geometry.integrationElement(r->position());
		for (int i = 0; i < size; i++)
			f[i] += fval[i] * (typename Function::ftype) weight * (typename Function::ftype) detjac;
#ifdef DEBUG
		std::cerr << "r =" << r->position() << std::endl;
		std::cerr << "weight =" << weight << " detjac =" << detjac << std::endl;
		std::cerr << "fval =" << fval << std::endl;
		std::cerr << "integral = " << f << std::endl;
#endif			

	}

	return f;
}

#endif
