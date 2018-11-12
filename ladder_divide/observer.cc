#include "itensor/all.h"
#include "hubbard_h/hubbard_d4_divide.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <complex>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <omp.h>
using namespace itensor;
#define PI 3.141592653589793
double get_wall_time()
{
    struct timeval time;
    if(gettimeofday(&time,NULL)){
        return 0;
    }

    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

int correlation_function(IQMPS psi, SiteSet sites, const std::string& StrOp1, const std::string& StrOp2, CMatrix& CorreFun, int Down, int begin, int end, const int N, const std::string& outFile)
{
//#pragma omp parallel for 
    int K;
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
        if(i >= Down+(begin*2-1) and i <= Down+(end*2-1))
        {
            K = N;
        }
        else if(i > N-4)
        {
            K = i+2;
        }
        else
        {
            K = i+4;
        }
        for(int j=i+2;j<=K;j+=2)
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
    N = 2*N; // two components
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
    auto ReadPsi = input.getYesNo("ReadPsi",false);
    auto ReadNum = input.getInt("ReadNum");
    auto WriteNum = input.getInt("WriteNum");
    
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    
    auto begin = input.getInt("begin");
    auto end = input.getInt("end");
    auto BudagBuPath = input.getString("BudagBuPath");
    auto BuBudagPath = input.getString("BuBudagPath");
    auto BddagBdPath = input.getString("BddagBdPath");
    auto BdBddagPath = input.getString("BdBddagPath");
    auto NupNupPath = input.getString("NupNupPath");
    auto SvNPath = input.getString("SvNPath");
    auto NupPath = input.getString("NupPath");
    auto NdnPath = input.getString("NdnPath");
    
    auto PBC = input.getYesNo("PBC");


    std::complex<double>im(0,1);
    println(sweeps);
    
    
    // Initialize the site degrees of freedom.
    //
    auto sites = HubbardD4Divide(N);
    readFromFile(format("sites_%d",ReadNum),sites);
        
    //
    IQMPS psi(sites);
    readFromFile(format("psi_%d",ReadNum),psi);
    Print(totalQN(psi));

//-------------  Correlation Functions ------------------
    double startTime = get_wall_time();
    CMatrix CorreFunBudagBu(N/2,N/2);
    auto StrOp1 = "Bupdag";
    auto StrOp2 = "Bup";
    auto Down = 0;
//omp_set_num_threads(4);
    println("Calculating BupdagBup");
    correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBudagBu, Down, begin, end, N, BudagBuPath);
    println("================================");
    double endTime = get_wall_time();
    println("BupdagBup Time : ", (double)(endTime-startTime), " s");
    

    startTime = get_wall_time();
    println("Calculating BupBupdag");
    CMatrix CorreFunBuBudag(N/2,N/2);
    StrOp1 = "Bup";
    StrOp2 = "Bupdag";
    correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBuBudag, Down, begin, end, N, BuBudagPath);
    println("================================");
    endTime = get_wall_time();
    println("BupBupdag Time : ", (double)(endTime-startTime), " s");
    
    startTime = get_wall_time();
    println("Calculating NupNup");
    CMatrix CorreFunNupNup(N/2,N/2);
    StrOp1 = "Nup";
    StrOp2 = "Nup";
    correlation_function(psi, sites, StrOp1, StrOp2, CorreFunNupNup, Down, begin, end, N, NupNupPath);
    println("================================");
    endTime = get_wall_time();
    println("NupNup Time : ", (double)(endTime-startTime), " s");
    //CMatrix CorreFunBddagBd(N/2,N/2);
    //StrOp1 = "Bdndag";
    //StrOp2 = "Bdn";
    //Down = 1;
    //correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBddagBd, Down, begin, end, N, BddagBdPath);
    
    //CMatrix CorreFunBdBddag(N/2,N/2);
    //StrOp1 = "Bdn";
    //StrOp2 = "Bdndag";
    //correlation_function(psi, sites, StrOp1, StrOp2, CorreFunBdBddag, Down, begin, end, N, BdBddagPath);
    


    // Print the final energy reported by DMRG
    return 0;
    }
