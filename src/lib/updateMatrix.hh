#ifndef DUNE_UPDATE_MATRIX_HH
# define DUNE_UPDATE_MATRIX_HH

#include <dune/common/dynmatrix.hh>

//#define DEBUG
#include "integrateEntityDynamic.hh"
#include "pointSource.hh"

//#define DEBUG

// A non-zero matrix element which have been change by the source
class sourceMatrixElement {

	public:
		sourceMatrixElement (long long unsigned int n, long long unsigned int m, double val) : n_(n), m_(m), value(val) {};
		long long unsigned int n () const {return n_;};
		long long unsigned int m () const {return m_;};
		double val () const {return value;};

	private:
		long long unsigned int n_,m_;
		double value;
};

template<class GridView, typename ftype, int fdim, class Source> 
std::vector<class sourceMatrixElement*>
updateMatrix (Dune::DynamicMatrix<ftype> & M, const GridView & gv, std::vector<Source*> slist)
{

	std::vector<class sourceMatrixElement*> smel;

	if (M.size() != gv.size (GridView::dimension) * (unsigned int) fdim)
	{
		std::cerr << "updateMatrix(): error: cannot update, matrix size wrong" << std::endl;
		return smel;
	}
#ifdef DEBUG	
	std::cerr << "updateMatrix(): M = " << M << std::endl;
#endif	

//#define DO_CHECKS
	for (typename std::vector<Source*>::iterator ps = slist.begin(); ps != slist.end(); ++ps)
	{
		Source * is = *ps;
		Dune::DynamicVector<ftype> m =  (*is) ();
		int vertexsize = is->entityVertexSize ();

#ifdef DO_CHECKS
		unsigned int size = vertexsize * fdim;
		size *= size;

		if (m.size() != size)
			DUNE_THROW (Dune::Exception, "Point Sources: User provided function returns DynamicVector of wrong size");
#endif
		for (int i = 0; i < vertexsize; i++)
			for (int j = 0; j < vertexsize; j++)
				for (int a = 0; a < fdim; a++)
					for (int b = 0; b < fdim; b++)
					{
						int ind = i * fdim * fdim * vertexsize + j * fdim * fdim + a * fdim + b;

						int ind1 = is->entityIndex (i) * fdim  + a;
						int ind2 = is->entityIndex(j) * fdim  + b;

						M[ind1][ind2] += m[ind];
						if (m[ind] != 0.0)
						{
#ifdef DEBUG						
							std::cerr << " updating M[" << ind1 << "][" << ind2 << "] by m[" 
								<< a << "][" << b << "]=" << m[ind] << std::endl;
#endif							
							smel.push_back (new sourceMatrixElement(ind1, ind2, m[ind]));
						}
#ifdef DEBUG					
						std::cerr << " += " << m[ind] 
							<< " = " << M[ind1][ind2] << std::endl;
#endif							
					}
	}

#ifdef DEBUG					
	std::cerr << "updateMatrix(): M = " << M << std::endl;
#endif

	return smel;
}
#undef DO_CHECKS

template<class GridView, typename ftype, int fdim, class Source> 
Dune::DynamicMatrix<ftype> assembleMatrixSources (const GridView & gv, std::vector<Source*> slist)
{
	Dune::DynamicMatrix<ftype> M;
	int N = gv.size (GridView::dimension) * fdim;
	M.resize (N,N);
	M = 0;

	updateMatrix<GridView, ftype, fdim, Source> (M, gv, slist);

	return M;
}

// Subtract the matrix elements due to the sources (previously saved in smel)
template<typename ftype> 
void 
revertMatrix (Dune::DynamicMatrix<ftype> & M, std::vector<sourceMatrixElement*> smel)
{

	for (typename std::vector<sourceMatrixElement*>::iterator me = smel.begin(); me != smel.end(); ++me)
	{
		long long unsigned int ind1 = (*me)->n();
		long long unsigned int ind2 = (*me)->m();

		assert (ind1 < M.mat_cols());
		assert (ind1 < M.mat_rows());

		M[ind1][ind2] -= (*me)->val();;
#ifdef DEBUG					
		std::cerr << " subtracting M[" << ind1 << "][" << ind2 << "] -= " << (*me)->val() << std::endl;
#endif							
	}

#ifdef DEBUG					
	std::cerr << "revertMatrix(): M = " << M << std::endl;
#endif

	return;
}
#undef DEBUG
#undef DO_CHECKS

#endif

