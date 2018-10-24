//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARDD3Divide_H
#define __ITENSOR_HUBBARDD3Divide_H
#include "itensor/mps/siteset.h"
#define sqrt2 1.414213562373
namespace itensor {

class HubbardSiteD3Divide;

using HubbardD3Divide = BasicSiteSet<HubbardSiteD3Divide>;

class HubbardSiteD3Divide
    {
    IQIndex s;
    public:

    HubbardSiteD3Divide() { }

    HubbardSiteD3Divide(IQIndex I) : s(I) { }

    HubbardSiteD3Divide(int n, Args const& args = Args::global())
        {
        auto conserveNf = args.getBool("ConserveNf",true);
        auto conserveSz = args.getBool("ConserveSz",true);
        int Up = (conserveSz ? +1 : 0),
            Dn = -Up;
        if(n%2 == 1) // up
        {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Nb=",0),
                    Index(nameint("Up ",n),1,Site),  QN("Sz=",1,"Nb=",1),
                    Index(nameint("UU ",n),1,Site), QN("Sz=", 2,"Nb=",2)};

        }
        else
        {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Nb=",0),
                    Index(nameint("Dn ",n),1,Site), QN("Sz=", -1,"Nb=",1),
                    Index(nameint("DD ",n),1,Site),  QN("Sz=",-2,"Nb=",2)};

        }
        /*
        if(conserveNf)
            {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Nf=",0),
                    Index(nameint("Up ",n),1,Site),  QN("Sz=",Up,"Nf=",1),
                    Index(nameint("Dn ",n),1,Site),  QN("Sz=",Dn,"Nf=",1),
                    Index(nameint("UpDn ",n),1,Site),QN("Sz=", 0,"Nf=",2)};
            }
        else //don't conserve Nf, only fermion parity
            {
            if(!conserveSz) Error("One of ConserveSz or ConserveNf must be true for Hubbard sites");

             s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Pf=",0),
                    Index(nameint("Up ",n),1,Site),  QN("Sz=",+1,"Pf=",1),
                    Index(nameint("Dn ",n),1,Site),  QN("Sz=",-1,"Pf=",1),
                    Index(nameint("UpDn ",n),1,Site),QN("Sz=", 0,"Pf=",0)};
            }
        */
        }

    IQIndex
    index() const { return s; }

    IQIndexVal
    state(std::string const& state)
        {
        if(state == "Emp") 
            {
            return s(1);
            }
        else 
        if(state == "Dn" || state == "Up") 
            {
            return s(2);
            }
        else 
        if(state == "DD" || state == "UU") 
            {
            return s(3);
            }
        else
            {
            Error("State " + state + " not recognized");
            }
        return IQIndexVal{};
        }

    IQTensor
    op(std::string const& opname,
       Args const& args) const
        {
        auto sP = prime(s);

        IQIndexVal Em(s(1)),
                   EmP(sP(1)),
                   Up(s(2)),
                   UpP(sP(2)),
                   Dn(s(2)),
                   DnP(sP(2)),
                   UU(s(3)),
                   UUP(sP(3)),
                   DD(s(3)),
                   DDP(sP(3));

        IQTensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up,UpP,1);
            Op.set(UU,UUP,2);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn,DnP,1);
            Op.set(DD,DDP,2);
            }
        else
        if(opname == "Bup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Up,UUP,sqrt2); 
            }
        else
        if(opname == "Bupdag")
            {
            Op.set(Up,EmP,1); 
            Op.set(UU,UpP,sqrt2); 
            }
        else
        if(opname == "Bdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Dn,DDP,sqrt2); 
            }
        else
        if(opname == "Bdndag")
            {
            Op.set(Dn,EmP,1); 
            Op.set(DD,DnP,sqrt2); 
            }
        else
            {
            Error("Operator \"" + opname + "\" name not recognized");
            }

        return Op;
        }
    };


} //namespace itensor

#endif
