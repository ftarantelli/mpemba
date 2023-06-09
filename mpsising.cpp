#include "mpsising.hpp"
#include <fstream>




int main(int argc, char *argv[])
    {
    
    auto expcutoff(12.);
    
    // The system will be a 2x20 ladder
    int Nx = 8;
    int Ny = 1;
    
    auto t = 0.05;
    auto tend = 1.;
    
    auto t0 = 0.001;
    
    auto J = 4.;
    auto g = 1.*J/2.;
    auto h = 0.5*J/2., hf = 0.*J/2.;
    auto sigma(10.), sigmaf(0.);

    std::string lineterm = "";
	for (int i=1;i<argc;i++)
		lineterm.append(std::string(argv[i]).append(" "));
    
while( argc > 1 ) {

	switch(argv[1][0]) {
		case 'N':
				Nx = atoi( &argv[1][1] );
			break;
		case 'g':
				g = atof( &argv[1][1] )*J/2.;
			break;
		case 'h':
			if(argv[1][1] == 'i')		h = atof( &argv[1][2] )*J/2.;
			else if(argv[1][1] == 'f')	hf = atof( &argv[1][2] )*J/2.;
			else						h = atof( &argv[1][1] )*J/2.;
			break;
		case 't':
			if(argv[1][1] == 'o')		t0 = atof( &argv[1][2] );
			else if(argv[1][1] == 'd')	tend = atof( &argv[1][2] );
			else						t = atof( &argv[1][1] );			
			break;
		case 'c':
				expcutoff = -atof( &argv[1][1] );
			break;
        case 's':
				if(argv[1][1] == 'R')	sigma = -std::sqrt(atof( &argv[1][2] ));
				else if(argv[1][1] == 'f')	sigmaf = atof( &argv[1][2] );
				else			sigma = atof( &argv[1][1] );
			break;
		case '-':
			//	std::sprintf(folder, "%s", &argv[1][1]);
			break;
		default:
			std::cerr << "Unlucky: Retry input values\n";
			exit (8);
	}
	++argv;
	--argc;
}
    tend = tend * double(Nx);
    int nsw = tend/t;
    int N = Nx*Ny;

    char output[30];

	std::sprintf(output, "data/mpsI%.3fF%.0fl%d.dat",  sigma, sigmaf, N );
	std::ofstream out_file(output, std::ios::out | std::ios::trunc);
    out_file.precision(16);

	//###chrono###chrono###chrono###chrono###chrono###
	char runlog[15];
	std::sprintf(runlog, "run.log");
	std::ofstream run_log(runlog, std::ios::out | std::ios::app);
	auto start = std::chrono::system_clock::now();
	std::time_t start_time = std::chrono::system_clock::to_time_t(start);
	run_log << output << " with INPUT: " << lineterm << "	at	" << std::ctime(&start_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###

    h = sigma / std::pow(N, 15./8.)*J/2.;
	hf = sigmaf / std::pow(N, 15./8.)*J/2.;

    // Make N spin 1/2's
    auto sites = SpinHalf(N, {"ConserveQNs=", false});

    // Make the Hamiltonian for rung-decoupled Heisenberg ladder
    auto ampo = AutoMPO(sites);    
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += -J, "Sx",i,"Sx",i+1;
        ampo += -g, "Sz",i;
        ampo += -h, "Sx",i;
        }
        ampo += -g, "Sz",N;
        ampo += -h, "Sx",N;
        
    auto H = toMPO(ampo);
    //printfln("Maximum bond dimension of H is %d",maxLinkDim(H));
    
    
    auto sweeps0 = Sweeps(5); //number of sweeps is 5
    sweeps0.maxdim() = 20,50,200,200,400;
    sweeps0.cutoff() = 1E-12;

    auto psi0 = randomMPS(sites);

    auto [energy0,psi] = dmrg(H,psi0,sweeps0,{"Silent=",true});
    //std::cout << energy0 << "\n";
    auto wfs = std::vector<MPS>(1);
    wfs.at(0) = psi0;

    //auto [en1,psi1] = dmrg(H,wfs,randomMPS(sites),sweeps0,{"Quiet=",true,"Weight=",20.0});

    //auto DeltaL = en1 - energy0;
    auto obsite = 3;
    
    auto psi1 = MPS(psi);
    auto psi2 = psi1;
    
    //out_file << 0. << "	" << real(measure(1, psi, sites)) << "	" << obsite << "\n";
	//out_file << 0. << "	" << real(measure(2, psi, sites)) << "	" << obsite << "\n";
    auto obs(0.);
	for(int nsite=1; nsite<=N; ++nsite)
            obs += real(measure(nsite, psi2, sites));
    out_file << 0. << "	" << obs/double(N) << "\n";
	//out_file << 0. << "	" << real(measure(4, psi, sites)) << "	" << obsite << "\n";
    // start TDVP, either one site or two site algorithm can be used by adjusting the "NumCenter" argument
    //println("----------------------------------------GSE-TDVP---------------------------------------");

    auto energy = real(innerC(psi1,H,psi1));
    printfln("Initial energy = %.5f", energy);
    
    
    
    // QUENCH
    h = h - hf;
    for(int i = 1; i <= N-1; ++ i)
        {
        ampo += h, "Sx",i;
        }
        ampo += h, "Sx",N;
    

    auto sweeps = Sweeps(1);
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = 1E-12;
    /*
    sweeps.maxdim() = 2000;
    sweeps.cutoff() = std::pow(1, -expcutoff);//1E-12;
    sweeps.niter() = 10;
    */

//std::cout << "hghgh11g\n";
    /*
    for(int n = 1; n <= nsw; ++n)
        {
        if(n < 3)
            {
            // Global subspace expansion
            std::vector<Real> epsilonK = {1E-12, 1E-12};
            addBasis(psi1,H,epsilonK,{"Cutoff",1E-8,
                                      "Method","DensityMatrix",
                                      "KrylovOrd",3,
                                      "DoNormalize",true,
                                      "Quiet",true});
            }
        
        // TDVP sweep
        energy = tdvp(psi1,H,-t,sweeps,{"Truncate",true,
                                        "DoNormalize",true,
                                        "Quiet",true,
                                        "NumCenter",1});
        }
    */
    //printfln("\nEnergy after imaginary time evolution = %.10f",energy);
    //printfln("Using overlap = %.10f", real(innerC(psi1,H,psi1)) );

    //println("-------------------------------------MPO W^I 2nd order---------------------------------------");
//std::cout << "hghghg\n";
    auto expH1 = toExpH(ampo,-(1-1_i)/2*t0*Cplx_i);
    auto expH2 = toExpH(ampo,-(1+1_i)/2*t0*Cplx_i);
    //printfln("Maximum bond dimension of expH1 is %d",maxLinkDim(expH1));
    auto args = Args("Method=","DensityMatrix","Cutoff=",1E-12,"MaxDim=",2000);
    out_file << "#   " << lineterm << '\n' << "# t				\n" ;

    for(int n = 1; n <= int(nsw*std::real(t/t0)); ++n)
        {
        psi2 = applyMPO(expH1,psi2,args);
        psi2.noPrime();
        psi2 = applyMPO(expH2,psi2,args);
        psi2.noPrime().normalize();

        if(n%int(std::real(t/t0)) == 0)
            {
            //printfln("\nMaximum bond dimension at time %.1f is %d ", n*t0, maxLinkDim(psi2));
            //printfln("Energy using overlap at time %.1f is %.10f", n*t0, real(innerC(psi2,H,psi2)) );
            obs = 0.;
            for(int nsite=1; nsite<=N; ++nsite)
                obs += real(measure(nsite, psi2, sites));
            out_file << n*t0 << "	" << obs/double(N) << "\n" << std::flush;
            }


        }

    //###chrono###chrono###chrono###chrono###chrono###
	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);
	run_log << output << " END with " << lineterm << "   " << std::ctime(&end_time) << std::flush;
	//###chrono###chrono###chrono###chrono###chrono###

    return 0;
    }
