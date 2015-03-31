//#ifdef HAVE_CONFIG_H
//#include "config.h"     
//#endif

//#define GRIDTYPE ALBERTAGRID
//#define GRIDDIM 3
//#define ALBERTA_DIM GRIDDIM 
//#define ALBERTAGRID

#include "config.h"
//#include <dune/grid/utility/gridtype.hh>

#include <iostream>
#include <dune/common/parallel/mpihelper.hh> // An initializer of MPI
#include <dune/common/exceptions.hh> // We use exceptions

#include <dune/grid/albertagrid.hh> //  Alberta  Grid
#include <dune/grid/albertagrid/dgfparser.hh>

//#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
//#include <dune/grid/io/file/dgfparser/dgfparser.hh>

#include <dune/grid/common/gridinfo.hh>
#include <dune/common/fvector.hh> // Field Vector
#include <dune/common/dynmatrix.hh> // Field Vector

//#define DEBUG
#include "assembleMatrix.hh"
#include "assembleVector.hh"
#include "saveField.hh"
#include "saveField2.hh"

#include "laplacian.hh"
#include "gamma.hh"
#include "sourceFunctions.hh"
#include "integrateField.hh"

#include "evalField.hh"

#include "setDirichlet.hh"
#include "setRobin.hh"
#include "calcSurfArea.hh"

#include "boundary-robin.hh"
#include "boundary-dirichelt.hh"

#ifdef DEBUG
# undef DEBUG
#endif

template<class GridView, typename ft, int fdim>
class entityFunctionVector : public entityP1Function<GridView, ft, fdim>
{
	public:
		entityFunctionVector (int i) {}

		template<class Entity> int size (const Entity & e) const
		{
				return fdim * this->vertexSize (e);
		}

		template<class Entity> Dune::DynamicVector<ft> 
			operator() (const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& x, 
					const Dune::FieldVector<typename GridView::ctype,GridView::dimension>& X, // global
					const Entity & e) const
			{
				Dune::DynamicVector<ft> V ;
				//int vertexsize = this->vertexSize(e);
				V.resize (size(e));
				V = 0;
/*
				for (int i = 0; i < vertexsize; i++)
					for (int a = 0; a < fdim; a++)
					{
						int ind2 = i * fdim + a;
						//V[ind2] =  this->basis[i].evaluateFunction(x) * x;
						V[ind2] =  0.;

					}
*/
#ifdef DEBUG
				std::cerr << "v = " << V << std::endl;
#endif			
				return  V;
			}
	protected:
};

int main(int argc, char** argv)
{

	double ks1 = 10;
	double ks2 = 2.;

	double L = 4.; /* surface-to-surface distance between enzymes */
	double H = 50.; /* box size */
	double Ro = 3.; /* size of enzymes (the same for E1 and E2) */

	if (argc >= 2)
	{
		L = atof (*(argv+1));
		std::cerr << "L=" << L << std::endl;
	}

	if (argc >= 3)
	{
		ks1 = atof (*(argv+2));
		std::cerr << "ks=" << ks1 << std::endl;
	}
	if (argc >= 4)
	{
		ks2 = atof (*(argv+3));
		std::cerr << "ks2=" << ks2 << std::endl;
	}

	try{
		//Maybe initialize Mpi
		Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
		std::cerr << "Hello World! This is test." << std::endl;

		if(Dune::MPIHelper::isFake)
			std::cerr<< "This is a sequential program." << std::endl;
		else
			std::cerr<<"I am rank "<<helper.rank()<<" of "<<helper.size()
				<<" processes!"<<std::endl;

#if (ENABLE_ALBERTA) && ((ALBERTA_DIM==3))

		static const int dim = 3;
		//const char * gridfile = "3dgrid.al";

		//typedef Dune::GridSelector::GridType Grid;
		//Dune::GridPtr<Grid> gridPtr ( "mesh.dgf" );

		typedef Dune::AlbertaGrid<dim,dim> Grid;
		Dune::GridPtr<Grid> gridPtr ( "mesh.dgf" );

		Grid & grid = * gridPtr ;

		//typedef Dune::AlbertaGrid<dim,dim> Grid;
		//Grid grid (gridfile);

		//grid.globalRefine (7);
		//Dune::gridinfo (grid);

		typedef typename Grid::LeafGridView GridView;
		const GridView & gv = grid.leafView();

		// list of points to finer the mesh
		std::vector<Dune::FieldVector<typename GridView::ctype, GridView::dimension> > plist;

		// This is the `source' -- first enzyme -- need here for finering the mesh
		Dune::FieldVector<typename GridView::ctype, GridView::dimension> xE1;
		xE1[0] = 0.5 * (H - L) - Ro; xE1[1] = 15.; xE1[2] = 0.;
		plist.push_back (xE1);
		Dune::FieldVector<typename GridView::ctype, GridView::dimension> xE2;
		xE2[0] = 0.5 * (H + L) + Ro; xE2[1] = 15.; xE2[2] = 0.;
		plist.push_back (xE2);

		const int fdim = 2;

		//
		// Matrix 
		//

		// Laplacian contribution to Galerkin matrix
		typedef entityFunctionLaplacian<GridView, double, fdim>  MIntegrand;
		MIntegrand mfunc (1);

		// List of point sources / enzymes 
		// 
		// how fast they move
		double dt = 1;

		//std::vector<int> smlist; // FIXME: dummy type
		//Dune::DynamicMatrix<double> M = assembleMatrix<GridView, MIntegrand> (gv, mfunc, 10, smlist);
		Dune::DynamicMatrix<double> M = assembleMatrix<GridView, MIntegrand> (gv, mfunc, 10);

		//
		// Vector contribution
		//
		typedef entityFunctionVector<GridView, double, fdim>  VIntegrand;
		VIntegrand vfunc (1);
		Dune::DynamicVector<double> V = assembleVector<GridView, VIntegrand> (gv, vfunc, 10);
		//std::cout << "V = " << V << std::endl;

		//
		// Set Robin BC
		//
	
		typedef RobinBC<typename GridView::ctype, GridView::dimension - 1 , double, fdim> Robin;
				
		// first enzyme
		Dune::FieldVector<double, fdim> u (0);
		u[0] = ks1;
		std::cerr << "xE1 = " << xE1 << std::endl;
		Robin bc1 (xE1, Ro, u);
		setRobin<GridView, double, fdim, Robin> (gv, M, V, bc1);

	
		// second enzyme
		Dune::FieldMatrix<double, fdim, fdim> eta (0);
		eta[0][0] = ks2;
		eta[1][0] = - ks2;
		std::cerr << "xE2 = " << xE2 << std::endl;
		Robin bc2 (xE2, Ro, eta); 
		setRobin<GridView, double, fdim, Robin> (gv, M, V, bc2);

		
		// Check surface areas
		std::cerr << "Area of enzyme 1: " <<  calcSurfaceArea (gv, bc1) << std::endl;
		std::cerr << "Area of enzyme 2: " <<  calcSurfaceArea (gv, bc2) << std::endl;

		//
		// set Dirichlet on cube sides
		//
		typedef DirichletBC<typename GridView::ctype, GridView::dimension - 1 , double, fdim> Dirichlet;
		const Dirichlet & bc = Dirichlet (); 
		setDirichlet<GridView, double, fdim, Dirichlet> (gv, M, V, 0, bc);

		//
		// Solve stationary problem
		//
		// Dune::DynamicVector<double> X;
		// X.resize (V.size());
		// M.solve (X, V);

		std::stringstream skE1;
		skE1 << std::setprecision(5) << ks1;

		std::stringstream skE2;
		skE2 << std::setprecision(5) << ks2;

		std::stringstream sL;
		sL << std::setprecision(5) << L;

		std::string fname = std::string ("kE1_" + skE1.str() + "-kE2_" + skE2.str() + "-L" + sL.str());

		typedef entityFunctionGamma<GridView, double, fdim>  GIntegrand;
		GIntegrand gfunc (1);
		Dune::DynamicMatrix<double> G = assembleMatrix<GridView, GIntegrand> (gv, gfunc, 10);

		double t = 0., tmax = 1000.;

		// M -> B = G + dt * M
		M *= dt;
		M += G;
		//std::cout << "M="  << M << std::endl;

		// file name 
		std::ostream * stream = new std::ofstream (fname + ".dat");

		// solution (rho) at t
		Dune::DynamicVector<double> rt;
		rt.resize (V.size());
		rt = 0.; // initial solution

		Dune::DynamicVector<double> B;
		B.resize (V.size());
		//G.mv (rt, B); // rt is zero any ways

		V *= dt;
		B += V;
		//std::cout << "b="  << b << std::endl;

		//M.solve (rt, B);
		M.invert();
		//M.mv (B, rt);
		//t = dt;

		unsigned int skip = 10;
		unsigned int iskip = 0;

		// to control how much they eat 
		double Ms = 0.;
		// mass eaten by an enzyme
		double ME = 0.;
//		typedef evalField<GridView, double, fdim>  EvalField;
//		EvalField f (gv, rt);

		//saveFieldToFile<Grid, double, fdim> (grid, rt, fname + "_t0.nodes");
/*
		*stream << t << "     " << integrateField<Grid, double, fdim> (grid, rt, 6)
			<< "   " << Ms << "    " << ME 
			//<< "   " << system.getNParticles() 
			<< std::endl;
*/
		while (t < tmax)
		{
		
			std::stringstream ss;
			ss << std::setprecision(5) << t;

			if ( !(iskip < skip))
			{
				iskip = 0;
		//		std::cerr << " -> saving profile" << std::endl;
				std::cerr << "t=" << t << std::endl;

				//saveFieldToFile<Grid, double, fdim> (grid, rt, fname + "_" + std::to_string(tint) + ".dat");
				//saveFieldToFile<Grid, double, fdim> (grid, rt, fname + "_t" + ss.str() + ".nodes");
			
				Dune::FieldVector<double, dim> x0 (0.0);
				Dune::FieldVector<double, dim> L (50);
				Dune::FieldVector<int, dim> N (101);

				//saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, L, N, fname + "_t" + ss.str() + ".dat");

		
				//double minval = 0.25, maxval = 19.5;
				double minval = 0.25;

		/*
				// x = 0 
				L[0] = 0.; N[0] = 0; x0[0] = 0.; 
				x0[1] = minval; L[1] = maxval; N[1] = 51.;
				x0[2] = minval; L[2] = maxval; N[2] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, L, N, fname + "_t" + ss.str() + "_x0.dat");
		*/
				// z = 0 
				L[2] = 0.; N[2] = 0; x0[2] = 0; 
				x0[0] = minval; L[0] = 49.75; N[0] = 101.;
				x0[1] = minval; L[1] = 29.75; N[1] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, L, N, fname + "_t" + ss.str() + "_z0.dat");
		/*
				// y = 20 
				L[1] = 0.; N[1] = 0; x0[1] = 20.; 
				x0[0] = minval; L[0] = maxval; N[0] = 51.;
				x0[2] = minval; L[2] = maxval; N[2] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, L, N, fname + "_t" + ss.str() + "_y20.dat");
		*/	
				//saveFieldToFile<Grid, double, fdim> (grid, rt, fname + ".nodes");
				//std::cerr << "rt=" << rt << std::endl;
			}
			iskip++;
			
			/*
			ME = .0;
			for (int i = 0; i < system.getNParticles(); i++)
			{	
				(system.getParticle(i))->Eat(f(system.getParticle(i)->position()));
				ME += (system.getParticle(i))->Eaten();
			}
			*/

			// Mass created by the emitting surface (see setRobin.hh)
			//Ms += dt * 5. * 0.2;
			
			//*stream << t << "      " << integrateField<Grid, double, fdim> (grid, rt, 6) << std::endl;
			*stream << t << "      " << integrateField<Grid, double, fdim> (grid, rt, 6)
			//	<< "   " << Ms << "    " << ME 
				<< std::endl;

			t += dt;

			
			// solver code
			//G.mv (rt, B);
			//B += V;
			//M.solve (rt, B);
			//
			// more efficient multiply code
			G.mv (rt, B);
			B += V;
			M.mv (B, rt);

		}

		//
		// Free lists
		//
		//for (typename std::vector<PSMatrix*>::iterator ps = smlist.begin(); ps != smlist.end(); ++ps)
		//	delete *ps;

#else
		std::cerr << "No 3D Alberta grid" << std::endl;
#endif

		return 0;
	}

	catch (Dune::Exception &e){
		std::cerr << "Dune reported error: " << e << std::endl;
	}
	catch (...){
		std::cerr << "Unknown exception thrown!" << std::endl;
	}

}
