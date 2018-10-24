//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_HUBBARDD4Divide_H
#define __ITENSOR_HUBBARDD4Divide_H
#include "itensor/mps/siteset.h"
#define sqrt2 1.414213562373
#define sqrt3 1.732050807569
namespace itensor {

class HubbardSiteD4Divide;

using HubbardD4Divide = BasicSiteSet<HubbardSiteD4Divide>;

class HubbardSiteD4Divide
    {
    IQIndex s;
    public:

    HubbardSiteD4Divide() { }

    HubbardSiteD4Divide(IQIndex I) : s(I) { }

    HubbardSiteD4Divide(int n, Args const& args = Args::global())
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
                    Index(nameint("UU ",n),1,Site), QN("Sz=", 2,"Nb=",2),
                    Index(nameint("U3 ",n),1,Site), QN("Sz=", 3,"Nb=",3)};

        }
        else
        {
            s = IQIndex{nameint("site=",n),
                    Index(nameint("Emp ",n),1,Site), QN("Sz=", 0,"Nb=",0),
                    Index(nameint("Dn ",n),1,Site), QN("Sz=", -1,"Nb=",1),
                    Index(nameint("DD ",n),1,Site),  QN("Sz=",-2,"Nb=",2),
                    Index(nameint("D3 ",n),1,Site),  QN("Sz=",-3,"Nb=",3)};

        }
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
        if(state == "D3" || state == "U3") 
            {
            return s(4);
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
                   DDP(sP(3)),
                   U3(s(4)),
                   U3P(sP(4)),
                   D3(s(4)),
                   D3P(sP(4));

        IQTensor Op(dag(s),sP);

        if(opname == "Nup")
            {
            Op.set(Up,UpP,1);
            Op.set(UU,UUP,2);
            Op.set(U3,U3P,3);
            }
        else
        if(opname == "Ndn")
            {
            Op.set(Dn,DnP,1);
            Op.set(DD,DDP,2);
            Op.set(D3,D3P,3);
            }
        else
        if(opname == "Bup")
            {
            Op.set(Em,UpP,1); 
            Op.set(Up,UUP,sqrt2); 
            Op.set(UU,U3P,sqrt3); 
            }
        else
        if(opname == "Bupdag")
            {
            Op.set(Up,EmP,1); 
            Op.set(UU,UpP,sqrt2); 
            Op.set(U3,UUP,sqrt3); 
            }
        else
        if(opname == "Bdn")
            {
            Op.set(Em,DnP,1); 
            Op.set(Dn,DDP,sqrt2); 
            Op.set(DD,D3P,sqrt3); 
            }
        else
        if(opname == "Bdndag")
            {
            Op.set(Dn,EmP,1); 
            Op.set(DD,DnP,sqrt2); 
            Op.set(D3,DDP,sqrt3); 
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
