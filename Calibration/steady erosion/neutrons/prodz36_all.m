% Function for calculating low-energy neutron parameters for a given set of inputs.
% Based on prodz36 function from Marrero et al. (2016).

function out=prodz36_all(pp,sp36,cp,sf);
% Where:
% pp = physical parameters, sp36 = sample parameters, cp = computed parameters,
% sf = scaling factors

%%% First, set up some parameters %%%

%depth independent variables for muon-produced neutrons
Pmu0=cp.a; %muonogenic neutron production at surface
Rmu=Pmu0/(sf(5)*pp.Pf0*cp.Reth);
Rmup=Rmu*pp.pEtha/cp.pEthss;

% Get Ca and K surface production from muons
ProdmuCa=cp.muon36(1,1);
ProdmuK=cp.muon36(2,1);

% parameters
Lambdafe = cp.Lambdafe;
Lth = cp.Lthss;
Leth = cp.Lethss;
Lambdamu = cp.b;

% correct for radiogenic neutrons
N_cor = cp.N36m - cp.N36r;

% neutron capture production pre-factors
k1 =  sf(5)*cp.Phistarethss*(1-cp.pEthss)*cp.feth/cp.Lambdaethss;
k2 =  sf(5)*(1+Rmu.*cp.Reth).*cp.FDeltaPhistareth*(1-cp.pEthss)*cp.feth/cp.Lambdaethss;
k3 =  sf(5)*Rmu.*cp.Phistarethss*(1-cp.pEthss)*cp.feth/cp.Lambdaethss;
k4 =  sf(5)*cp.Phistarthss*cp.fth/cp.Lambdathss;
k5 =  sf(5)*(1+Rmup).*cp.SFDeltaPhistarethss*cp.fth/cp.Lambdathss;
k6 =  sf(5)*(1+Rmup.*cp.Rth).*cp.SFDeltaPhistarthss*cp.fth/cp.Lambdathss;
k7 =  sf(5)*Rmup.*cp.Phistarthss*cp.fth/cp.Lambdathss;

% muon production
k8 = ProdmuCa+ProdmuK;

% Other outputs
out.Lambdafe = Lambdafe;
out.Lambdamu = Lambdamu;
out.Lth = Lth;
out.Leth = Leth;
out.k1 = k1;
out.k2 = k2;
out.k3 = k3;
out.k4 = k4;
out.k5 = k5;
out.k6 = k6;
out.k7 = k7;
out.k8 = k8;
out.R = cp.N36r;


