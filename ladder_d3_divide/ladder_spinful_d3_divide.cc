#include "itensor/all.h"
#include "itensor/mps/sites/hubbard_d3_divide.h"
#include <cmath>
#include <complex>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
using namespace itensor;
#define PI 3.14159265357
double get_wall_time()
{
	struct timeval time;
	if(gettimeofday(&time,NULL)){
		return 0;
	}

	return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

int correlation_function(IQMPS psi, SiteSet sites, const std::string& StrOp1, const std::string& StrOp2, CMatrix& CorreFun, int Down, const int N, const std::string& outFile)
{
#pragma omp parallel for 
	for(int i=Down+1;i<=N-2;i+=2)
	{
		IQMPS newPsi = psi;
		newPsi.position(i);
		auto Op1 = sites.op(StrOp1,i); 
		IQTensor C = newPsi.A(i);
		C *= Op1;
		auto ir = commonIndex(newPsi.A(i),newPsi.A(i+1),Link);
		C *= dag(prime(newPsi.A(i),Site,ir));
		//C *= dag(prime(prime(psi.A(i),Site),ir)); // A*Op1*A_star 
		for(int j=i+2;j<=N;j+=2)
		{
			auto Op2 = sites.op(StrOp2,j);
			if(j == i+2)
			{
				C *= newPsi.A(j-1);
				C *= dag(prime(newPsi.A(j-1),Link));
				
				IQTensor CC = C*newPsi.A(j);
				CC *= Op2;
				auto jl = commonIndex(newPsi.A(j),newPsi.A(j-1),Link);
				CC *= dag(prime(newPsi.A(j),jl,Site));
				//CC *= dag(prime(prime(psi.A(j),Site),jl));
				CorreFun((i-1)/2,(j-1)/2) = CC.cplx();
			}
			else
			{
				C *= newPsi.A(j-2);
				C *= dag(prime(newPsi.A(j-2),Link));
				C *= newPsi.A(j-1);
				C *= dag(prime(newPsi.A(j-1),Link));
				
				ITensor CC = C*newPsi.A(j);
				CC *= Op2;
				auto il = commonIndex(newPsi.A(j),newPsi.A(j-1),Link);
				CC *= dag(prime(prime(newPsi.A(j),Site),il));
				CorreFun((i-1)/2,(j-1)/2) = CC.cplx();
			}

		}
	}
	std::ofstream outfile(outFile);
	for(int i=1;i<=N/2-1;i++)
	{
		for(int j=i+1;j<=N/2;j++)
		{
			outfile << i << "," << j << " " << CorreFun(i-1,j-1) << std::endl;
		}
	}

	return 0;
}

int main(int argc, char* argv[])
    {
//--------------- parameters --------------
    if(argc != 2) { printfln("Usage: %s inputfile",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");
    
    auto N = input.getInt("N");
    auto Nup = input.getInt("Nup");
    auto Ndn = input.getInt("Ndn");
    N = 2*N;
    auto nsweeps = input.getInt("nsweeps");
    auto Jx = input.getReal("Jx");
    auto Jy = input.getReal("Jy");
    auto Uuu = input.getReal("Uuu");
    auto Udd = input.getReal("Udd");
    auto Uud = input.getReal("Uud");
    auto phi = input.getReal("phi");
	phi = phi*PI;
    auto quiet = input.getYesNo("quiet",false);
    int writem = input.getYesNo("writem",false);
    
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
	
	auto BudagBuPath = input.getString("BudagBuPath");
	auto BuBudagPath = input.getString("BuBudagPath");
	auto BddagBdPath = input.getString("BddagBdPath");
	auto BdBddagPath = input.getString("BdBddagPath");
	auto SvNPath = input.getString("SvNPath");
	auto NupPath = input.getString("NupPath");
	auto NdnPath = input.getString("NdnPath");
	
	auto PBC = input.getYesNo("PBC");
	//auto Nhalf = 6;
	//auto N = 2*Nhalf;  
    //auto sweeps = Sweeps(20);
	//sweeps.maxm() = 1000;
	//sweeps.minm() = 100;
	//sweeps.cutoff() = 1E-14;
    //auto Jx = 1.0;
    //auto Jy = 1.0;
    //auto U = 2.0;
	//auto phi = 1.0/4.0*2*PI;
	//auto phi = 1.57;
	//auto quiet = true;


	std::complex<double>im(0,1);
    println(sweeps);
	
    
    // Initialize the site degrees of freedom.
    //
    auto sites = HubbardD3Divide(N);
    //
    // Create the Hamiltonian using AutoMPO
    //
    auto ampo = AutoMPO(sites);
    for(int i = 1; i <= N-7; i+=4) 
        {
		ampo += -Jx,"Bupdag",i,"Bup",i+4;
		ampo += -Jx,"Bupdag",i+4,"Bup",i;
		ampo += -Jx,"Bdndag",i+1,"Bdn",i+5;
		ampo += -Jx,"Bdndag",i+5,"Bdn",i+1;
		
		ampo += -Jx,"Bupdag",i+2,"Bup",i+6;
		ampo += -Jx,"Bupdag",i+6,"Bup",i+2;
		ampo += -Jx,"Bdndag",i+3,"Bdn",i+7;
		ampo += -Jx,"Bdndag",i+7,"Bdn",i+3;
        }
	if(PBC)
	{
		// PBC
		ampo += -Jx,"Bupdag",1,"Bup",N-3;
		ampo += -Jx,"Bupdag",N-3,"Bup",1;
		ampo += -Jx,"Bdndag",2,"Bdn",N-2;
		ampo += -Jx,"Bdndag",N-2,"Bdn",2;
   
		ampo += -Jx,"Bupdag",3,"Bup",N-1;
		ampo += -Jx,"Bupdag",N-1,"Bup",3;
		ampo += -Jx,"Bdndag",4,"Bdn",N;
		ampo += -Jx,"Bdndag",N,"Bdn",4;
	}
	for(int i = 1; i<=N-1; i+=2)
		{
        ampo += Uuu/2,"Nup",i,"Nup",i;
		ampo += -Uuu/2,"Nup",i;
        ampo += Udd/2,"Ndn",i+1,"Ndn",i+1;
		ampo += -Udd/2,"Ndn",i+1;
		ampo += -Uud/2,"Nup",i,"Ndn",i+1;
		}
    for(int i = 1; i <= N-3; i += 4)
        {
		auto r = (i+3)/4.0;
        ampo += -Jy*exp(im*r*phi),"Bupdag",i,"Bup",i+2;
        ampo += -Jy*exp(-im*r*phi),"Bupdag",i+2,"Bup",i;
        ampo += -Jy*exp(im*r*phi),"Bdndag",i+1,"Bdn",i+3;
        ampo += -Jy*exp(-im*r*phi),"Bdndag",i+3,"Bdn",i+1;
        }
    auto H = IQMPO(ampo);
    //
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    //
    auto state = InitState(sites);
    int p = Nup;
    int q = Ndn;
	
	for(int i = 1;i <= N/2-1;i+=2)
    {
    	if(p == 0 and q == 0)
    	{
    		state.set(i,"Emp");
    		state.set(i+1,"Emp");
    		state.set(N-i,"Emp");
    		state.set(N-i+1,"Emp");
    	}
		else
		if(p == 0)
		{
			state.set(i,"Emp");
			state.set(N-i,"Emp");
		}
		else
		if(q == 0)
		{
			state.set(i+1,"Emp");
			state.set(N-i+1,"Emp");
		}
    	if(p > 1)
    	{
    		println("Site ",i," Up");
    		println("Site ",N-i," Up");
    		state.set(i,"Up");
    		state.set(N-i,"Up");
    		p -= 2;
    	}
    	else
    	if(p > 0)
    	{
    		println("Site ",i," Up");
    		state.set(i,"Up");
    		state.set(N-i,"Emp");
    		p -= 1;
    	}
    	
    	if(q > 1)
    	{
    		println("Site ",i+1," Dn");
    		println("Site ",N-i+1," Dn");
    		state.set(i+1,"Dn");
    		state.set(N-i+1,"Dn");
    		q -= 2;
    	}
    	else
    	if(q > 0)
    	{
    		println("Site ",i+1," Dn");
    		state.set(i+1,"Dn");
    		state.set(N-i+1,"Emp");
    		q -= 1;
    	}
    
	}
    auto psi = IQMPS(state);

    Print(totalQN(psi));

    //
    // Begin the DMRG calculation
    //
    auto energy = dmrg(psi,H,sweeps,{"Quiet",quiet,"WriteM",writem});

    //
    // Measure densities
    //
    //Vector upd(N),dnd(N);
	auto halfN = N/2;
    Vector densityUp(halfN);
    Vector densityDn(halfN);
	auto k = 0;
    for(int j = 1; j <= N-1; j += 2)
        {
        psi.position(j);
		densityUp(k) = (dag(prime(psi.A(j),Site))*sites.op("Nup",j)*psi.A(j)).real();
		psi.position(j+1);
		densityDn(k) = (dag(prime(psi.A(j+1),Site))*sites.op("Ndn",j+1)*psi.A(j+1)).real();
        k++;
		}
	

	std::ofstream outfile(NupPath);
    println("Density Nup:");
    for(int j = 0; j < halfN; ++j)
	{
        printfln("%d %.10f",1+j,densityUp(j));
		outfile << j+1 << " " << densityUp(j) << std::endl;
	}
	outfile.close();
    println();
	
	std::ofstream outfileNdn(NdnPath);
    println("Density Ndn:");
    for(int j = 0; j < halfN; ++j)
	{
        printfln("%d %.10f",1+j,densityDn(j));
		outfileNdn << j+1 << " " << densityDn(j) << std::endl;
	}
	outfileNdn.close();
    println();

	
//------------- von Neumann entanglement ------------------
	Vector SvN(N/2-1);
	for(int i=2;i<N;i+=2)
	{
		psi.position(i);
		IQTensor twf = psi.A(i)*psi.A(i+1);
		auto U = psi.A(i);
		IQTensor S,V;
		auto spectrum = svd(twf,U,S,V);
		Real Sv = 0;
		for(auto p : spectrum.eigs())
		{
			if(p > 1E-12) Sv += -p*log(p);
		}
		SvN((i-1)/2) = Sv;
	}

	println("von Neumann entanglement entropy");
	for(int i=1;i<N/2;i++)
	{
		println("Bond",i," ",SvN(i-1));
	}
	
	std::ofstream outfileSvN(SvNPath);
	for(int i=1;i<N/2;i++)
	{
		outfileSvN << i << " " << SvN(i-1) << std::endl;
	}
	outfileSvN.close();

//-------------  Correlation Functions ------------------
	double startTime = get_wall_time();
	omp_set_num_threads(4);
	CMatrix CorreFunBudagBu(N/2,N/2);
	auto StrOp1 = "Bupdag";
	auto StrOp2 = "Bup";
	auto Down = 0;
	
	correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBudagBu, Down, N, BudagBuPath);
	
	
	CMatrix CorreFunBuBudag(N/2,N/2);
	StrOp1 = "Bup";
	StrOp2 = "Bupdag";
	correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBuBudag, Down, N, BuBudagPath);
	
	CMatrix CorreFunBddagBd(N/2,N/2);
	StrOp1 = "Bdndag";
	StrOp2 = "Bdn";
	Down = 1;
	correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBddagBd, Down, N, BddagBdPath);
	
	CMatrix CorreFunBdBddag(N/2,N/2);
	StrOp1 = "Bdn";
	StrOp2 = "Bdndag";
	correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBdBddag, Down, N, BdBddagPath);
	
	double endTime = get_wall_time();
	println("Correlation Functions Time : ", (double)(endTime-startTime), " s");
	println("Correlation Functions:");
	
	//for(int i=1;i<=N/2-1;i++)
	//{
	//	for(int j=i+1;j<=N/2;j++)
	//	{
	//		println("(",i,",",j,")",CorreFunBudagBu(i-1,j-1));
	//		println();
	//	}
	//}



    // Print the final energy reported by DMRG
    //
    printfln("\nGround State Energy = %.10f",energy);
    return 0;
    }
