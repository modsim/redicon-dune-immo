//#ifdef HAVE_CONFIG_H
//#include "config.h"     
//#endif

//#define GRIDTYPE ALBERTAGRID
//#define GRIDDIM 3
//#define ALBERTA_DIM GRIDDIM 
//#define ALBERTAGRID

#include "config.h"
//#include <dune/grid/utility/gridtype.hh>
//#define _GNU_SOURCE
//#ifdef HAVE_GETOPT_H
#  include <getopt.h>
//#endif /* HAVE_GETOPT_H */

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

#include "gamma.hh"
#include "sourceFunctions.hh"
#include "integrateField.hh"

#include "evalField.hh"

#include "setDirichlet.hh"
#include "setRobin.hh"
#include "calcSurfArea.hh"

#include "bridge-entityFunction.hh"
#include "bridge-robinBC.hh"
#include "bridge-forces.hh"

#include "boundary-dirichelt.hh"

#ifdef DEBUG
# undef DEBUG
#endif

const char * myname = "immpoplate-bridge";
#define MAINTAINER   "valiska@gmail.com"

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

// Diagonal (isotropic) diffusion matrix for fdim fields in dim dimension
template<typename ct, int dim, typename ft, int fdim>
class Diffusion
{
	public:
		Diffusion (const Dune::FieldVector<ft, fdim> & D) : D(D) {};

		Dune::FieldMatrix<ft, dim * fdim, dim * fdim> 
			operator() (const Dune::FieldVector<ct,dim>& x) const
			{
				Dune::FieldMatrix<ft, dim * fdim, dim * fdim> m (0.0);
				
				for (int a = 0; a < fdim; a++)
					for (int i = 0; i < dim; i++)
						m[a*dim + i][a*dim + i] = D[a];

				return  m;
			};
	protected:
		const Dune::FieldVector<ft, fdim> & D;
};


int main(int argc, char** argv)
{
	/* this to say nmol/s ? */
	//double ks1 = 1.64;
	//double ks2 = .04;
	// this is as calculated
	//double ks1 = 0.000003;
	//double ks2 = 0.00008;
	// just play:
	double ks1 = 0.3;
	double ks2 = 8.;

	double L = 4.; /* surface-to-surface distance between enzymes */
	double H[3] = {50., 50., 30.}; /* box size */

	double Ro1 = 5.; /* size of the 1st enzymes */
	double Ro2 = 2.5; /* size of the 2nd enzymes */
	double bridge_width = 2. * Ro2;
	double bridge_amplitude = 0.047;
	std::string type;
	std::string fname;

	/* 
	 * Command-line parser 
	 */  
	int c =0;

#define ARGS    "hvL:A:W:1:2:H:T:o:"
	while (c != EOF) 
	{
#ifdef HAVE_GETOPT_LONG
		int option_index = 0;
		static struct option long_options[] = 
		{
			{"help", no_argument, NULL, 'h'},
			{"version", no_argument, NULL, 'v'},

			{"enzyme-separation", required_argument, NULL, 'L'},
			{"bridge-amplitude", required_argument, NULL, 'A'},
			{"bridge-width", required_argument, NULL, 'W'},
	//		{"bridge-type", required_argument, NULL, 'T'},

			{"box-size", required_argument, NULL, 'H'},
			{"output", required_argument, NULL, 'o'},

			{"enyzme-1-rate", no_argument, NULL, '1'},
			{"enyzme-2-rate", no_argument, NULL, '2'},

		};
		switch ((c = getopt_long (argc, argv, ARGS, long_options, &option_index))) 
#else 
			switch ((c = getopt (argc, argv, ARGS))) 
#endif /* HAVE_GETOPT_LONG */
			{

				case '1': 
					ks1 = atof (optarg);
					std::cerr << "ks1=" << ks1 << std::endl;
					break;

				case '2': 
					ks2 = atof (optarg);
					std::cerr << "ks2=" << ks1 << std::endl;
					break;

				case 'L': 
					L = atof (optarg);
					std::cerr << "L=" << L << std::endl;
					break;

				case 'T': 
					type = std::string (optarg);
					break;

				case 'o': 
					fname = std::string (optarg);
					break;

				case 'H':
				{
					std::stringstream ss(optarg);
					int i = 0;
					double val = 0.;
					while (ss >> val)
					{
						if (i < 3)
							H[i] = val;
						else
						{
							std::cerr << "Error: Use immoplata-bridge -H Hx,Hy,Hz" << std::endl;
							return 0;
						}
						if (ss.peek() == ',')
							ss.ignore();
					}
				}
					break;

				case 'W': 
					bridge_width = atof (optarg);
					std::cerr << "W=" << bridge_width << std::endl;
					break;

				case 'A': 
					bridge_amplitude = atof (optarg);
					std::cerr << "A=" << bridge_amplitude << std::endl;
					break;

				case 'h': /* help */
					fprintf (stderr,
							"Usage: %s [options]\n"
							"Where options are:\n"
							"  -H, --box=Hx,Hy,Hz                       box size (default: %g,%g,%g).\n"

							"  -L, --enzyme-separation=VAL              distance between enzymes (default: %g).\n"
							"  -A, --bridge-amplitude=VAL               amplitude of the bridge force (default: %g).\n"
							"  -W, --bridge-width=VAL                   width of the bridge (default: %g).\n"
							"  -1, --enzyme-1-rate=VAL                  rate of enzyme 1 producing intermediate (default: %g).\n"
							"  -2, --enzyme-2-rate=VAL                  rate of enzyme 2 transforming intermediate (default: %g).\n"
							"  -o, --output=FILE                        output file name.\n"
							"  -v, --version                            display version information and exit.\n"
							"  -h, --help                               display this help and exit.\n"
							"Report bugs to %s\n", myname, H[0], H[1], H[2], L, bridge_amplitude, bridge_width, ks1, ks2, MAINTAINER);
					return 0; /* success */
					break;

				case '?': /* wrong options */
					fprintf (stderr, "Try `%s --help' for more information.\n", myname);
					return 1; /* failure */
			}
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

		typedef Dune::AlbertaGrid<dim,dim> Grid;
		Dune::GridPtr<Grid> gridPtr ( "mesh.dgf" );

		Grid & grid = * gridPtr ;

		typedef typename Grid::LeafGridView GridView;
		const GridView & gv = grid.leafView();

		const int fdim = 2;

		//
		// Matrix 
		//

		// Diffusion
		const Dune::FieldVector<double, fdim> D (1.);
		typedef Diffusion<GridView::ctype, GridView::dimension, double, fdim>  Diff;
		Diff diff (D);

		//
		// Force
		//
		// force between the plate and particles
		//typedef ForceElectrostatic<GridView::ctype, GridView::dimension, double, fdim>  Force;
		typedef ForceLJ<GridView::ctype, GridView::dimension, double, fdim>  Force;

		// charge product (q_bridge * q_intermediate)
		Dune::FieldVector<double, fdim> Q (bridge_amplitude);
		Q[1] = 0.0; // products do not interact with the bridge
		Force force (Q, .5); // external force (non integrated, parcile-particle)

		// plate geometry:
		// 	position
		Dune::FieldVector<double, 3> x0;
		for (int i = 0; i < 2; i++)
			x0[i] = 0.5 * H[i]; 
		// x axis
		x0[0] += 0.5 * (Ro2 - Ro1);
		// z axis
		x0[2] = 0.;

		// half size
		Dune::FieldVector<double, 2> w;
		w[0] =  0.5 * L + Ro1 + Ro2;
		w[1] = 0.5 * bridge_width;

		//typedef RectBridge<double, 3, Force> BridgeForce;
		//BridgeForce F (x0, w, force, 60); // force from the the whole 'bridge' acting on intermediate molecules

		//
		// This is Morse force, we need no integraton
		//
		typedef MorseBridge<double, 3, double, 2> BridgeForce;
		BridgeForce F (x0, w, bridge_amplitude, 0.05); // force from the the whole 'bridge' acting on intermediate molecules

		//
		// Current in the Galerkin matrix (Force)
		//
		typedef entityFunctionForce<GridView, double, fdim, Diff, BridgeForce>  MIntegrand;
		MIntegrand mfunc (D,F);

		double dt = 0.1;

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
	
		typedef RobinBCForce<typename GridView::ctype, GridView::dimension - 1 , double, fdim, BridgeForce> Robin;
				
		// Enzymes

		// first enzyme
		// position
		Dune::FieldVector<typename GridView::ctype, GridView::dimension> xE1;
		xE1[0] = 0.5 * (H[0] - L) - Ro1;
		xE1[1] = 0.5 * H[1]; 
		xE1[2] = 0.;
		// strength
		Dune::FieldVector<double, fdim> u (0);
		u[0] = ks1;
		std::cerr << "xE1 = " << xE1 << std::endl;
		Robin bc1 (xE1, Ro1, u, F);
		setRobin<GridView, double, fdim, Robin> (gv, M, V, bc1);

	
		// second enzyme
		// strength
		Dune::FieldMatrix<double, fdim, fdim> eta (0);
		eta[0][0] = ks2;
		eta[1][0] = - ks2;
		// position
		Dune::FieldVector<typename GridView::ctype, GridView::dimension> xE2;
		xE2[0] = 0.5 * (H[0] + L) + Ro2; 
		xE2[1] = 0.5 * H[1]; 
		xE2[2] = 0.;
		std::cerr << "xE2 = " << xE2 << std::endl;
		Robin bc2 (xE2, Ro2, eta, F); 
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

		if (fname.empty())
		{
			std::stringstream skE1;
			skE1 << std::setprecision(5) << ks1;

			std::stringstream skE2;
			skE2 << std::setprecision(5) << ks2;

			std::stringstream sL;
			sL << std::setprecision(5) << L;

			fname = std::string ("kE1_" + skE1.str() + "-kE2_" + skE2.str() + "-L" + sL.str());
		}

		typedef entityFunctionGamma<GridView, double, fdim>  GIntegrand;
		GIntegrand gfunc (1);
		Dune::DynamicMatrix<double> G = assembleMatrix<GridView, GIntegrand> (gv, gfunc, 10);

		double t = 0., tmax = 100.;

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
		//double Ms = 0.;
		// mass eaten by an enzyme
		//double ME = 0.;
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
			
				Dune::FieldVector<double, dim> x0 (0.0); // starting point (from:)
				Dune::FieldVector<double, dim> Size (50); // size (to:)
				Dune::FieldVector<int, dim> N (101); // number of points
			
				//saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, L, N, fname + "_t" + ss.str() + ".dat");

		
				//double minval = 0.25, maxval = 19.5;
				//double minval = 0.25;

				double eps = 0.25;
				// x = 0 
				Size[0] = 0.; N[0] = 0; x0[0] = 0.; 
				x0[1] = eps; Size[1] = H[1] - eps; N[1] = 51.;
				x0[2] = eps; Size[2] = H[2] - eps; N[2] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, Size, N, fname + "_t" + ss.str() + "_x0.dat");
		
				// z = 0 
				Size[2] = 0.; N[2] = 0; x0[2] = 0; 
				x0[0] = eps; Size[0] = H[0] - eps; N[0] = 51.;
				x0[1] = eps; Size[1] = H[1] - eps; N[1] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, Size, N, fname + "_t" + ss.str() + "_z0.dat");
		
				// y = max
				Size[1] = 0.; N[1] = 0; x0[1] = H[1]; 
				std::stringstream Hss;
				Hss << std::setprecision(5) << x0[1];

				x0[0] = eps; Size[0] = H[0] - eps; N[0] = 51.;
				x0[2] = eps; Size[2] = H[2] - eps; N[2] = 51.;
				saveFieldToFile2<Grid, double, fdim> (grid, rt, x0, Size, N, fname + "_t" + ss.str() + "_y" + Hss.str() + ".dat");
			
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

		delete stream;

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
