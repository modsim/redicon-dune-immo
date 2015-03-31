#ifndef DUNE_POINT_SOURCE_HH
#define DUNE_POINT_SOURCE_HH

#include "dune/common/dynvector.hh" 
#include "diracDelta.hh"

template<class GridView, class Function> 
class PointSource 
{
	public:
		enum {dimension = GridView::dimension};
		typedef typename GridView::ctype ctype;
		typedef typename GridView ::template Codim<0>::Entity Entity;

		PointSource (const GridView & gv, diracDelta<ctype, dimension> & delta, const Function & f, int key)
			: gv(gv), delta (delta), function(f), key(key), entityptr (nullptr)
		{
			if (!findEntity())
				DUNE_THROW (Dune::Exception, "Cannot determine my entity in grid, probably not in grid");
		}

		PointSource (const GridView & gv, const Dune::FieldVector<ctype, dimension> & pos, const Function & f, int key)
			: gv(gv), delta (pos), function(f), key(key), entityptr (nullptr)
		{
			if (!findEntity())
				DUNE_THROW (Dune::Exception, "Cannot determine my entity in grid, probably not in grid");
		}

		~PointSource () {}

		bool findEntity()
		{
			for (typename GridView :: template Codim<0> :: Iterator it = gv. template begin <0>(); it != gv. template end<0>(); ++it)
				if (delta.inside (it))
				{
					entityptr = &(*it);
					local_ = it->geometry().local(delta.x());
					break;
				}

			if (!entityptr)
				return false;

			const typename GridView::IndexSet & set = gv.indexSet();
			indices.resize (entityVertexSize());
			for (int i = 0; i < entityVertexSize(); i++)
				indices[i] = set.subIndex (*entityptr, i, dimension);

			return true;
		}

		bool resetPosition (const Dune::FieldVector<ctype, dimension> & y)
		{
			const Dune::FieldVector<ctype, dimension> & x_ = x();
			delta.resetPosition(y);
			if (findEntity())
			{
				//std::cerr << "position reset to " << y << std::endl;
				//std::cerr << "new incedes:";
				//for (int i = 0; i < entityVertexSize(); i++)
				//	std::cerr << " " << indices[i];
				//std::cerr << '\n';

				return true;
			}
			else
			{
				delta.resetPosition(x_);
				std::cerr << "resetPosition(): error: cannot reset as I cannot find an element" << std::endl;
				return false;
			}
		}

		const int getKey () const {return key;};

		const Dune::FieldVector<typename GridView::ctype, GridView::dimension> x () const {return this->delta.x();};
		const Dune::FieldVector<typename GridView::ctype, GridView::dimension> global () const {return x();};

		const Dune::FieldVector<typename GridView::ctype, GridView::dimension> local () const {return local_;};

		int entityIndex (int ivertex)
		{
			if (ivertex > entityVertexSize())
				DUNE_THROW (Dune::Exception, "vertex index out of range");
			return indices[ivertex];
		}


		// FIXME: is this a bug? shall we store the vertex size instead?
		int entityVertexSize () const
		{
			return function.vertexSize (*entityptr);
		}

		Dune::DynamicVector<typename Function::ftype> eval (void) const {return function (local(), *entityptr); }

		Dune::DynamicVector<typename Function::ftype> operator() (void) const
		{
	//		return function (x(), *entityptr); 
			// std::cerr << "f(x=" << local() <<  ")=" ; std::cerr << function (local(), *entityptr) << std::endl;
			return function (local(), *entityptr, key); 
		}

		Dune::DynamicVector<typename Function::ftype> 
		operator() (const Dune::FieldVector<typename GridView::ctype, GridView::dimension> & x) const
		{
			if (x == delta.x())
				return function (x, *entityptr, key); 
			else
			{
				Dune::DynamicVector<typename Function::ftype> r (Function::fdimension);
				r = 0.;
				return r;
			}
		}

	private:
		const GridView & gv;
		diracDelta<ctype, dimension> delta;
		const Function & function;
		const int key;
		Entity * entityptr; // FIXME: need this?
		Dune::FieldVector<ctype, dimension> local_; // local coordinates within enityt (element)
		std::vector<int> indices; 
};

#endif
