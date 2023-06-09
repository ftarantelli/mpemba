#include <cstring>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <math.h>

#include <armadillo>
// hi A&M
using namespace arma;
typedef std::complex<double> complex;

double g, hi, hf, h, tau;


class chain {
	public:
		int L, states;
		int PBC;
		Mat<int> spin;
		cx_dvec gs;
		sp_cx_dmat Ham, Ham0, Vp;
		double energy;
		double Mval, Aval, Zval;
		double time;
		
		chain(int _L, int _pbc) {
			L = _L;
			states = std::pow(2, _L);
			PBC = _pbc;
			spin.zeros(states, _L);
			gs.zeros(states);
			Ham0.set_size(states, states);
			Ham.set_size(states, states);
			Vp.set_size(states, states);
			energy = 0.;
			time = 0.;
			Mval = 0.;
			Zval = 0.;
			for(int i = 1; i <= states; ++i) {
				int z = i-1;
				for(int j=0; j < L; ++j) {
					spin(i-1, j) = int(std::pow(-1, z%2+1));
					z = int(z/2);
				}
			
			}
		}
		
		void BuildHam(void);
		void GroundState(void);
		void Measure(int point);
		
		void derivate( double t, const cx_dvec& _gs, cx_dvec& k);
		void RungeKuttaPsi(double _dt);
		
		void SuzukiTrotter(double _dt);
		
};


void chain::BuildHam(void) {

		
		double temp;
		for(int i = 1; i <= states; ++i) {
			temp = 0.;
			for(int j = 1; j < L; ++j) {
				int est1 = i - spin(i-1,j-1)*std::pow(2, j-1);
				//Ham(est1-1,i-1) = - h;
				Vp(est1-1,i-1) = - 1.;
				int est2 = est1-spin(est1-1,j+1-1)*std::pow(2, j);
				Ham0(est2-1,i-1) = - 1.;
				temp += - g * spin(i-1,j-1);
			}
			int est1 = i - spin(i-1,L-1)*std::pow(2, L-1);
			//Ham(est1-1,i-1) = - h;
			Vp(est1-1,i-1) = - 1.;
			temp += - g*spin(i-1,L-1);
			Ham0(i-1,i-1) = temp;
		}
		if(PBC==1) {
			for(int i = 1; i <= states; ++i) {
				int est1 = i - spin(i-1,0);
				int est2 = est1-spin(est1-1,L-1)*std::pow(2, L-1);
				Ham0(est2-1,i-1) = - 1.;
			}
		}
		//Ham = Ham0 + h*Vp;
}


void chain::GroundState(void) {

		//sp_dmat HH = real(Ham);
		//eigs_opts opts;
		//opts.maxiter = 1000;
		cx_dvec eigval;
		cx_dmat eigvec;

		//vec = eigs_sym(HH);
		eigs_gen(eigval, eigvec, Ham, 1, "sr");
		for(int i=0; i<states; ++i)
			gs(i) = (eigvec.col(0))(i);// complex( (eigvec.col(0))(i), 0.);
		gs = gs / norm(gs);

		//std::cout  << eigvec.size() <<"\n" ;

}



void chain::Measure(int point) {

	complex obs(0.), obs1(0.);
	assert(point > 0 && point <= L);

	for(int j=1; j <= L; ++j) {
		for(int i=1; i <= states; ++i) {
		
//std::cout << j << " " << i <<  " aa\n";
			//obs +=(gs(i-1)*gs.t()(i-1))*complex(spin(i-1, point-1),0.);
			int est = i - spin(i-1, j-1)*std::pow(2, (j-1));
			//std::cout << est << " " << gs <<  " aa\n";
			obs += conj(gs(i-1))*gs(est-1);
			obs1 += complex(spin(i-1, j-1))* conj(gs(i-1)) * gs(i-1);
		}
	}

	Mval = obs.real() / double(L);
	Zval = obs1.real() / double(L);
	/*
	for(int i=0; i < gs.size(); ++i)
		obs += ( gs(i) * gs.t()(i) ) * complex(spin(i, point-1), 0.);
	Mval = obs.real();
	*/
	//delete &obs;
}



void chain::derivate( double t, const cx_dvec& _gs, cx_dvec& k) {
	//int dim = _gs.size();
	//SparseMatrix<complex> Ham(states, states);
	//std::cout << _Ham  << "\n";
	//std::cout << "a01\n";
	//SparseMatrix<complex> HH(Ham.rows(), Ham.cols());
	//for(int i=0; i<Ham.cols(); ++i)
	//	for(int j=0; j<Ham.rows(); ++j)
	//		HH(i, j) = complex(Ham(i,j), 0.);
	k = complex(0., -1.) * ( Ham * _gs );
	//std::cout << "a02\n";
}

void chain::RungeKuttaPsi(double _dt) {

	int dim(gs.size());
	cx_dvec k1(dim), k2(dim), k3(dim), k4(dim), temp;
	double time(0.);
//std::cout << "aa\n";
	derivate( time         , gs            , k1);
//	std::cout << "aa1\n";
	temp = gs + k1*_dt/2.;
	derivate( time + _dt/2., temp, k2);
	temp = gs + k2*_dt/2.;
	derivate( time + _dt/2., temp, k3);
	temp = gs + k3*_dt;
	derivate( time + _dt   , temp, k4);

	gs = gs + _dt/6.*( k1 + 2.*k2 + 2.*k3 + k4 );
	//time = time + _dt;
	//gs.normalise();
	gs = gs / norm(gs);
	//delete &k1, &k2, &k3, &k3, &temp;

}





int main(int argc, char *argv[]) {

	int Lsize(12); int pbc(0);
	double dt(0.02), tmax(10.*double(Lsize));
	g = 1.; hi = -0.25; hf = 0.25;// tau = 100;
	double sigma = 5., sigmaf = 0.;	// , omega = 1.

	std::string lineterm = "";
	for (int i=1;i<argc;i++)
		lineterm.append(std::string(argv[i]).append(" "));

while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'L':
				Lsize = atoi( &argv[1][1] );
			break;
		case 'g':
				g = atof( &argv[1][1] );
			break;
		case 'h':
			if(argv[1][1] == 'i')		hi = atof( &argv[1][2] );
			else if(argv[1][1] == 'f')	hf = atof( &argv[1][2] );
			//else {	h = atof( &argv[1][1] );	h1 = h;  }
			break;
		case 'b':
				pbc = atoi( &argv[1][1] );
			break;
		case 'd':
				dt = atof( &argv[1][1] );
			break;
		case 't':
				if(argv[1][1] == 'l')	tmax = double(Lsize)*atof( &argv[1][2] );
				else tmax = atof( &argv[1][1] );
			break;
		//case 'o':
				//omega = atof( &argv[1][1] );
			//break;
		case 's':
				if(argv[1][1] == 'R')	sigma = -std::sqrt(atof( &argv[1][2] ));
				else if(argv[1][1] == 'f')	sigmaf = atof( &argv[1][2] );
				else			sigma = atof( &argv[1][1] );
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}


char output[30];

	std::sprintf(output, "data/ismpI%.3fF%.0fl%d.dat",  sigma, sigmaf, Lsize );
	std::ofstream out_file(output, std::ios::out | std::ios::trunc);
//	std::sprintf(output, "data/armkzUp%.0fSig%.0fL%dC%dA.dat", upsilon*10., 10*sigma, Lsize, n_cycle );
//	std::ofstream out_file1(output, std::ios::out | std::ios::trunc);
//}
//std::cout << log(datum::e) << "\n";

out_file.precision(16);


	//###chrono###chrono###chrono###chrono###chrono###
	char runlog[15];
	std::sprintf(runlog, "run.log");
	std::ofstream run_log(runlog, std::ios::out | std::ios::app);
	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	run_log << output << " with INPUT: " << lineterm << "	at	" << std::ctime(&start_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###

	out_file << "#   " << lineterm << '\n' << "# t/L		Mx		Mz		A		\n" ;

class chain C(Lsize, pbc);
class chain E(Lsize, pbc);

C.BuildHam();
//C.GroundState();
E = C;
cx_dvec temp_gs(C.states);

int events = int(tmax/dt);
std::cout << tmax << "  " << hf << "\n";
//exit(8);


	hi = sigma / std::pow(Lsize, 15./8.);
	hf = sigmaf / std::pow(Lsize, 15./8.);

	h = hi;
	//h = 0.5;
	C.Ham = C.Ham0 + hi * C.Vp;
	C.GroundState();
	temp_gs = C.gs;

	double Mtemp;

	//hi = std::abs(hi);
	C.Ham = C.Ham0 + hf * C.Vp;


for(int inx=0; inx <= events; ++inx) {
	if ( inx%(events/100) == 0) {
		//std::cout << int(inx/float(events)*100.) << "%\n";
		out_file << std::flush;
	}
	
		C.time += dt;
		//h = C.time / tau;
		//if ( int( (C.time-t00)/tau/2./hi ) % 2 == 0 )
		//if (sum + dir*dt/tau <= hf && dir==1)
		//	dir = 1.;
		//else
		//	dir = -1.;
		//sum += dir * dt / tau;
//std::cout << dir << " " << sum << "\n";
		//C.Ham = C.Ham + C.Vp * dir * (dt / tau);
		//h = 0.;

		//C.SuzukiTrotter(dt);
		C.RungeKuttaPsi(dt);
		//C.BuildHam(true);

		//E.time += dt;

		//if( sum + dir*dt/tau > hf ) {
		//if((sum + dt/tau) * std::pow(upsilon, 15./23.) * std::pow(Lsize, 15./8.) >= abs(sigma)) {
		C.Measure(3);

		if( abs(C.Mval-Mtemp) < std::pow(10., -7.)  ){
			std::cout << inx << "\n";
			break;
		}
		E.Ham = C.Ham;
		E.GroundState();
		out_file << C.time/double(Lsize) << "		" << C.Mval << "		" << C.Zval << "		" << norm((E.gs.t()*C.gs)) << '\n';// << std::flush;
		//}
		Mtemp = C.Mval;
	}

	//###chrono###chrono###chrono###chrono###chrono###
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###

/*
std::sprintf(output, "data/gsL%dg%.0lfhi%.0f_%.0lfhf%.0f_%.0lf.dat", Lsize, g, hi, 10.*abs(hi), hf, 10.*hf );
std::ofstream out_file1(output, std::ios::out | std::ios::trunc);
out_file << "#	time		Mval	\n";


	h = hi;


while( C.time  <  tmax ) {

	int frac=5;
	if( frac*(C.time-t00)/(tmax-t00) - int(frac*(C.time-t00)/(tmax-t00)) < frac*dt/(tmax-t00) )
		std::cout << int( (C.time - t00)/(tmax - t00) * 100 ) << "%\n";

	C.time += dt;
	h = C.time / tau;
	
	C.BuildHam();
	C.GroundState();

	C.Measure(3);	
	out_file1 << C.time << "		" << C.Mval << '\n';

}
*/

return(0.);
}























