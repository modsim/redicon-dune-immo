#ifndef DUNE_INSIDE_ENTITY_HH
# define DUNE_INSIDE_ENTITY_HH

#include <dune/common/fvector.hh>

//#define DEBUG
#define EPS_ENTITY 1.e-04

template<class Entity, typename ct, int dim> 
bool insideEntity (const Dune::FieldVector<ct, dim> & x, const Entity & e)
{ 
	//			const typename Entity::Geometry geometry = e.geometry();
	//			const Dune::GeometryType & gt = geometry.type();
	const Dune::GeometryType & gt = e->geometry().type();

	if ( gt.isCube() && (dim == 1) )
	{
		Dune::FieldVector<ct, dim> x1 =  e->geometry().corner(0);
		Dune::FieldVector<ct, dim> x2 =  e->geometry().corner(1);

		if ( (x[0] >= x1[0]) && (x[0] <= x2[0]) )
		{
			return true;
		}
	} 
	else if ( gt.isCube() && (dim == 2) )
	{
#ifdef DEBUG		
		for (int i = 0; i < e->geometry().corners(); i++)
			std::cerr << "corner " <<  i << ": " << e->geometry().corner(i) << std::endl;
		std::cerr << std::endl;
#endif
		// So it seems so:
		//  2 - 3
		//  |   |
		//  0 - 1
		double xmin = (e->geometry().corner(0))[0];
		double ymin = (e->geometry().corner(0))[1];
		double xmax = (e->geometry().corner(3))[0];
		double ymax = (e->geometry().corner(3))[1];
		if ( (x[0] >= xmin) && (x[0] <= xmax) 
				&& (x[1] >= ymin) && (x[1] <= ymax) )
		{
			//local_ = e->geometry().local(global_);
			return true;
		}
		else
			return false;
	}
	// this should be any dimension
	else if ( gt.isSimplex() )
	{
#ifdef DEBUG			
		std::cerr << "x= " <<  x << std::endl;
		for (int i = 0; i < e->geometry().corners(); i++)
			std::cerr << "corner " <<  i << ": " << e->geometry().corner(i) << std::endl;
		std::cerr << std::endl;
#endif

		// FIXME: rewrite this?
		int nvertices = e->geometry().corners();
#ifdef DEBUG		
		std::cerr << "nvertices: " << nvertices << std::endl;
#endif		
		if (nvertices != dim + 1)
		{
			std::cerr << "Number of vertices in a simplex confront with the dimensionality "
				<< nvertices << " vs " << dim << std::endl;
			return false;
		}
		Dune::FieldVector<ct, dim> Vi [nvertices];
		for (int i = 0; i < nvertices; i++)
		{
			Vi[i] = (e->geometry().corner(i));
#ifdef DEBUG			
			std::cerr << "vertex " <<  i << ": " << Vi[i] << std::endl;
#endif			
		}

		// fill in matrices and calculate determinants
		Dune::FieldMatrix <ct, dim+1, dim+1> D [nvertices+1];
		ct dets [nvertices+1];
		for (int m = 0; m < nvertices+1; m++)
		{
			for (int i = 0; i < nvertices; i++)
			{
				if ( i == m - 1 )
					for (int j = 0; j < dim; j++)
						(D[m])[i][j] = (x)[j];
				else
					for (int j = 0; j < dim; j++)
						(D[m])[i][j] = (Vi[i])[j];
				(D[m])[i][nvertices-1] = 1.;
			}
#ifdef DEBUG			
			std::cerr << "D" <<  m << ": " << D[m] << std::endl;
#endif
			dets[m] = (D[m]).determinant();
			// set to zero if it's zero, otherwise can lead to no inside below
			if ( fabs(dets[m]) < EPS_ENTITY )
				dets[m] = 0.0;
		}


		if ( fabs(dets[0]) < EPS_ENTITY)
		{
			std::cerr << "Simplex with zero determinant" << std::endl;
			return false;
		}

		ct detsum = 0.;
		bool inside = true;
		for (int m = 1; m < nvertices+1; m++)
		{
			detsum += dets[m];
			if ( (dets[0] * dets[m]) < 0.0 )
			{
#ifdef DEBUG			
				std::cerr << "*** Not inside! *** " << std::endl;
#endif				
				inside = false;
			}
		}

		if ( fabs(dets[0] - detsum) > EPS_ENTITY )
		{
			std::cerr << "Determinand check failed D0 = "
				<< dets[0] 
				<< " != \\sum_i Di = " 
				<< detsum
				<< std::endl;
			return false;
		}

#ifdef DEBUG
		if (inside)
			std::cerr << "Inside!" << std::endl;
#endif
#undef EPS_ENTITY

		return inside;
		//local_ = e->geometry().local(global_);
	}
	else
	{
		std::cerr << "geometry type: " << e->geometry().type() << std::endl;
		for (int i = 0; i < e->geometry().corners(); i++)
			std::cerr << "corner " <<  i << ": " << e->geometry().corner(i) << std::endl;
		std::cerr << std::endl;

		DUNE_THROW (Dune::Exception, "Geometry/dimension type not implemented");
	}

	return false;
};

//#undef DEBUG
#endif

