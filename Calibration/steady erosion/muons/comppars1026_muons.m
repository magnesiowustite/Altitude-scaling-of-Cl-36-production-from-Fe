% Code for calculating muon production parameters for Be-10
% modified from Marrero et al. (2016)
% - Muonfluxsato changed to 3-exponential term muon model and

function cp=comppars1026_muons(pp,sp36) %changed from sp1026 to sp36 11.8.22
%
% First, set maxdepth to a default value if not supplied.
%
cp.depthvector = [0,logspace(0,3.5,25)];
%
% Setup Lambdafe.
%
cp.Lambdafe=sp36.Lambdafe;
%
%---------------------------------------------------------------------
%
% The following are modifications for Sato/Heisinger muons.
%
%
% Pmunf
%
% old formulation:
%  cp.Pmunf=0.0000058*pp.Yf*pp.Phimuf0*exp(-cp.Zs/pp.Lambdamu);
%
%  New Formulation uses Greg's Heisinger code to calculate the fluxes at
% a vector of depths to create a table that can be referred to later

% RcEst = 14.9.*((cos(abs(d2r(sp36.latitude)))).^4); %Dipolar estimate for these purposes
% The RcEst above was originally used and has been replaced by Nat Lifton's
% new model.  He  fit the function below to trajectory-traced cutoffs for the
% modern dipolar field strength, so it should be more accurate overall.

dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];
%   Trajectory-traced dipolar estimate for these purposes
   RcEst = (dd(1)*cos(deg2rad(sp36.latitude)) + ...
       dd(2)*(cos(deg2rad(sp36.latitude))).^2 + ...
       dd(3)*(cos(deg2rad(sp36.latitude))).^3 + ...
       dd(4)*(cos(deg2rad(sp36.latitude))).^4 + ...
       dd(5)*(cos(deg2rad(sp36.latitude))).^5 + ...
       dd(6)*(cos(deg2rad(sp36.latitude))).^6);

  %set the constants to Be-10
  Natoms = pp.Natoms10;
  %  sigma190 = pp.sigma190_10;
  %    pp.delsigma190 = pp.delsigma190_10; % not used
  sigma0=pp.sigma010;
  k_negpartial = pp.k_negpartial10;
  fstar=pp.fstar10;
  %    pp.delk_neg = pp.delk_neg10; % not used

  flux=Muon_model_muons(cp.depthvector,sp36.P,RcEst,pp.SPhiInf);
%store the output fluxes that we need
%save coefficients and attenuation lengths
cp.a = flux.a;
cp.b = flux.b;
cp.c = flux.c;
cp.d = flux.d;
cp.e = flux.e;
cp.f = flux.f;

%parameters
  aalpha = 1.0;
  Beta = 1.0;

  Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*0)) + 50.7.*(1-exp(-5.05e-7.*0));

  % fast muon production of Be-10

  cp.P_fast10 = cp.e.*Beta.*(Ebar.^aalpha).*sigma0.*Natoms;

  % negative muon capture
  cp.P_neg10_a = cp.a.*k_negpartial.*fstar;
  cp.P_neg10_c = cp.c.*k_negpartial.*fstar;
  cp.P_negtotal10 = cp.P_neg10_a+cp.P_neg10_c;



