% Code for calculating production parameters modified from Marrero et al. (2016)

% Modified to accomodate analytical erosion rate calculations
% - Muonfluxsato changed to exponential muon model and
% - flux depth vector changed to flux0 for muon production on Ca and K.

% Generates a complete set of parameters for chlorine-36 dating
% from the basic physics parameters, sample specific parameters,
% and scale factors.

function cp=comppars36_NC(pp,sp36)

%Where:
% pp = physical parameters, sp36 = sample parameters 36

% Setup depthvector from 0 to 10,000 g/cm^2
cp.depthvector = [0,logspace(0,4,25)];

% Setup Lambdafe.
%
cp.Lambdafe=sp36.Lambdafe;
cp.ls=sp36.ls;

% Parameters for neutron-capture production, unmodified from original code

% Ni are the atoms per gram for C, Na, Mg, Al, Si, P, K, Ca, Ti,
% Mn, Fe, followed by Cl, B, Sm, Gd, U, Th, Li, Cr.  Note that C starts in
% entry 3 of the array to leave room for H(2) and O(1).
% Target elements will be calculated in cp.tni in a similar fashion.
% Target elements in order are K, Ca, Ti, Fe, Cl.

cp.ni=zeros(21,1);
for i=3:13
  cp.ni(i)=sp36.ci(i)*pp.table(i,9)*1.0e19;
end
for i=14:21
  cp.ni(i)=sp36.ci(i)*pp.table(i,9)*1.0e15;
end

cp.ni(2)=sp36.qtotal*0.111*pp.ag;
cp.ni(1)=cp.ni(2:13)'*pp.table(2:13,10);

cp.tni=zeros(5,1);
for i=1:3
  cp.tni(i)=sp36.tci(i)*pp.table((i+8),9)*1.0e19;
end
cp.tni(4)=sp36.tci(4)*pp.table(13,9)*1.0e19;
cp.tni(5)=sp36.tci(5)*pp.table(14,9)*1.0e15;
%
% Fraction of O, H, C, Na, Mg, etc.
%
cp.fi=zeros(21,1);
cp.fi=cp.ni .* pp.table(:,1)/ ...
    pp.ag;

% Ass is the average atomic weight of the subsurface material.
%
cp.Ass=pp.ag/sum(cp.ni);
%
% Sigma_sc,ss=sumproduct(Ni,ssc_i)/1.0e24
% Eqn 3.22? in Gosse & Phillips
%
cp.Sigmascss=cp.ni'*pp.table(:,3)/ ...
    1.0e24;
%
% Sigma_th,ss=sumproduct(Ni,ssth_i)/1.0e24
% Eqn 3.6 in Gosse & Phillips
%

cp.Sigmathss=cp.ni'*pp.table(:,4)/ ...
    1.0e24;
%
% I_eff=sumproduct(Ni,Ia_i)/1.0e24;
% I is the effective resonance integral for
% absorption of epithermal neutrons
% Eqn 3.9 in G&P
%
cp.Ieff=cp.ni'*pp.table(:,5)/1.0e24;
%
% X_ss=sumproduct(Ni,xi,ssc_i)/sumproduct(Ni,ssc_i)
%
cp.Xss=...
    sum((cp.ni.*pp.table(:,2)).*pp.table(:,3))/...
    (cp.ni'*pp.table(:,3));
%
% Sigma_eth,ss=Xss*(Ieff+Sigma_sc,ss)
%
cp.Sigmaethss=cp.Xss*...
    (cp.Ieff+cp.Sigmascss);
%
% Lambda_f,t=(4.3/3.3)*Lambda_f,a
%
%cp.Lambdaft=(4.3/3.3)*pp.Lambdafe;
%
% Lambda_eth,a=1/Sigma_eth,a
%
cp.Lambdaetha=1/pp.Sigmaetha;
%
% Lambda_eth,ss=1/Sigma_eth,ss
%
cp.Lambdaethss=1/cp.Sigmaethss;
%
% Lambda_th,ss=1/Sigma_th,ss
%
cp.Lambdathss=1/cp.Sigmathss;
%
% Added new variables for individual spallation rxns
%
cp.PsCa=cp.tni(2)*pp.PsCa0*(40.078/6.022E+23); % () is a conversion to at 36Cl/at Ca;
cp.PsK=cp.tni(1)*pp.PsK0*(39.0983/6.022E+23); % () is a conversion to at 36Cl/at K;
cp.PsTi=cp.tni(3)*pp.PsTi0;
cp.PsFe=cp.tni(4)*pp.PsFe0;

%
% B=sumproduct(xi,Ni,ssc_i)/1.0e24;
% D34 on Parameters page in CHLOE
% B is the scattering rate parameter
%
cp.B=...
    sum((cp.ni.*pp.table(:,2)).*pp.table(:,3))/1.0e24;
%
% p(Eth)_ss=exp(-Ieff/B)
%
cp.pEthss=exp(-cp.Ieff/cp.B);
%
% Reth=sqrt(Ass/14.5)
% Reth is the production rate of epithermal neutrons from
% fast neutrons in the air at land/atmosphere boundary
% P. 1498 in G&P
%
cp.Reth=sqrt(cp.Ass/14.5);
%
% Rth=pEthss/pEtha
% Eqn 3.20 in G&P
%
cp.Rth=cp.pEthss/pp.pEtha;
%
% D_eth,a=1/(3*Sigma_sc,a*(1-2/(3*Aa)))
% Eqn 3.16 in G&P
%
cp.Detha=1/(3*pp.Sigmasca*...
                        (1-2/(3*pp.Aa)));
cp.Da=cp.Detha;
%
% D_th,ss=(3*Sigma_sc,ss*(1-2/(3*Ass)))^(-1)
% Eqn 3.33 in G&P
%
cp.Dthss=(3*cp.Sigmascss*(1-2/(3*cp.Ass)))^(-1);
cp.Dethss=cp.Dthss;
%
% Lambda_th,a=1/Sigma_th,a
%
cp.Lambdatha=1/pp.Sigmatha;
%
% L_eth,a=1/sqrt(3*SSc_a*Seth_a)
% epithermal diffusion length (eqn 3.21? in G&P)
cp.Letha=1/sqrt(3*pp.Sigmasca* ...
                            pp.Sigmaetha);
%
% L_eth,ss=1/sqrt(Sigma_eth,ss*(3*Sigma_sc,ss))
%
cp.Lethss=1/sqrt(cp.Sigmaethss*...
                           (3*cp.Sigmascss));
%
% L_th,ss=sqrt(D_th,ss/Sigma_th,ss)
% thermal neutron diffusion length (eqn 3.34? in G&P)
%
cp.Lthss=sqrt(cp.Dthss/ ...
                          cp.Sigmathss);
%
% f_th=((s_th,cl*Ncl)/Sigma_th,ss)/1.0e24;
% Eqn 3.32 in G&P
%
cp.fth=((pp.table(14,4)*cp.tni(5))/...
    (cp.Sigmathss))/1.0e24;
%
% f_eth=((la_cl*Ncl)/Ieff)/1.0e24;
%
cp.feth=((pp.table(14,5)*cp.tni(5))/...
    cp.Ieff)/1.0e24;
%
% Phistaretha=Pf0/(Sigmaetha-Detha/Lambdafe^2)
% Eqn 3.26 in G&P
%
cp.Phistaretha=pp.Pf0/...
    (pp.Sigmaetha-cp.Detha/ ...
     cp.Lambdafe^2);
%
% Phistartha=p(Eth)_a*Sigma_eth,a*PhiStar_eth,a/(Sigma_th,a-Da/Lf_a^2)
% Eqn 3.38 in G&P
%
cp.Phistartha=pp.pEtha*pp.Sigmaetha*...
    cp.Phistaretha/...
    (pp.Sigmatha-cp.Da/cp.Lambdafe^2);
%
% Phistarethss=Pf0*Reth/(Sigmaethss-Dethss/Lambda_fa^2)
% Eqn 3.26 in G&P
%
cp.Phistarethss=pp.Pf0*cp.Reth/...
    (cp.Sigmaethss-cp.Dethss/ ...
     cp.Lambdafe^2);
%
% Phistarthss=pEtha*Sigmaethss*PhiStarethss*Rth/(Sigmathss-Dthss/Lambda_fa^2);
% Eqn 3.38 in G&P
%
cp.Phistarthss=pp.pEtha*cp.Sigmaethss*...
    cp.Phistarethss*cp.Rth/...
    (cp.Sigmathss-cp.Dthss/cp.Lambdafe^2);
%
% DeltaPhistareth=Phistaretha-Phistarethss
% Eqn 3.29 in G&P
%
cp.DeltaPhistareth=cp.Phistaretha-...
    cp.Phistarethss;
cp.DeltaPhistarethss=cp.DeltaPhistareth;
cp.DeltaPhistaretha=-cp.DeltaPhistareth;
%
% DeltaPhistarth=Phistartha-Phistarthss
% Eqn 3.3.41 in G&P
%
cp.DeltaPhistarth=cp.Phistartha-...
    cp.Phistarthss;
cp.DeltaPhistarthss=cp.DeltaPhistarth;
cp.DeltaPhistartha=-cp.DeltaPhistarth;
%
% DeltaPhistarstareth=Phistarethss-Detha*Phistaretha/Dethss
% Eqn 3.30 in G&P
%
cp.DeltaPhistarstareth=cp.Phistarethss-...
    cp.Detha*cp.Phistaretha/ ...
    cp.Dethss;
%
% FDeltaPhistarethss=((Detha*DelPhiStareth/Letha)-
%                     (Dethss*DelPhiDblStareth/Lfa))/
%                     ((Detha/Letha)+(Dethss/Lethss))
% Eqn 3.28 in G&P
%
cp.FDeltaPhistarethss=...
    (cp.Detha*cp.DeltaPhistareth/cp.Letha-...
     (cp.Dethss*cp.DeltaPhistarstareth/cp.Lambdafe))/...
    ((cp.Detha/cp.Letha)+...
     (cp.Dethss/cp.Lethss));
cp.FDeltaPhistareth=cp.FDeltaPhistarethss;
%
% FDeltaPhistaretha=((Dethss*(-DeltaPhistareth)/Lethss)-(Dethss*DeltaPhistarstareth/Lambdafe))/...
%  ((Detha/Letha)+Dethss/Lethss))
% Eqn 3.28 in G&P
%
cp.FDeltaPhistaretha=...
    ((cp.Dethss*(-cp.DeltaPhistareth)/cp.Lethss)-...
     (cp.Dethss*cp.DeltaPhistarstareth/cp.Lambdafe))/...
    ((cp.Detha/cp.Letha)+...
     (cp.Dethss/cp.Lethss));
%
% Eqn 3.39 in G&P
%
cp.SFDeltaPhistaretha=pp.pEtha*pp.Sigmaetha*cp.FDeltaPhistaretha/...
    (pp.Sigmatha-cp.Da/cp.Letha^2);
%
%Eqn 3.39 in G&P
%
cp.SFDeltaPhistarethss=cp.Sigmaethss*pp.pEtha*cp.FDeltaPhistareth*cp.Rth/...
    (cp.Sigmathss-cp.Dthss/cp.Lethss^2);
%
% Eqn 3.42 in G&P
%
cp.DeltaSFDeltaPhistareth=cp.SFDeltaPhistaretha- ...
    cp.SFDeltaPhistarethss;
cp.DeltaSFDeltaPhistaretha=-cp.DeltaSFDeltaPhistareth;
cp.DeltaSFDeltaPhistarethss=cp.DeltaSFDeltaPhistareth;
%
% cp.SFDeltaPhistartha.  Note correction to Eqn 3.40 in G&P
%
%
% Corrected version used below.
%
cp.SFDeltaPhistartha=...
    (cp.Da*(cp.Phistartha/cp.Lambdafe-cp.SFDeltaPhistaretha/cp.Letha)-...
     cp.Dthss*(cp.Phistarthss/cp.Lambdafe+cp.SFDeltaPhistarethss/cp.Lethss)+...
     (cp.Dthss/cp.Lthss)*(-cp.DeltaPhistarth-cp.DeltaSFDeltaPhistareth))/...
    (cp.Dthss/cp.Lthss+cp.Da/pp.La);
%
% SFDeltaPhistarthss.
%
%
% Corrected version of G&P eqn 3.40 used below.
%
cp.SFDeltaPhistarthss=...
    (cp.Da*(cp.Phistartha/cp.Lambdafe-cp.SFDeltaPhistaretha/cp.Letha)-...
     cp.Dthss*(cp.Phistarthss/cp.Lambdafe+cp.SFDeltaPhistarethss/cp.Lethss)+...
     (cp.Da/pp.La)*(cp.DeltaPhistarthss+cp.DeltaSFDeltaPhistarethss))/...
    (cp.Dthss/cp.Lthss+cp.Da/pp.La);
%
% Zs
%
%cp.Zs=sp36.ls*sp36.rb;
%
%---------------------------------------------------------------------
%
% The following are modifications for Sato/Heisinger muons.
% Modified 8.10.22 to parameterize muons for analytical calculations

%  New Formulation uses Greg's Heisinger code to calculate the fluxes at
% a vector of depths to create a table that can be referred to later
% Actually get the full data from P_mu_total -- needed later
    % Use middle of sample
%     mu = P_mu_total((sample.thick.*sample.rho./2),sample.pressure,mconsts,'yes');
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

   cp.RcEst = RcEst; %set as a variable to feed to muon_model
     %     % %Use mean Solar Modulation Parameter (SPhiInf)

%the following assignments of pp are so the code doesn't crash.  We don't
%use the output values, so they are meaningless. If they are included this
%way, we don't have to recode the muon code for each nuclide. Muons are
%time-independent as coded.

flux=Muon_model(cp.depthvector,sp36.P,RcEst,pp.SPhiInf);
%store the output fluxes that we need
%save pre-factor and attenuation length
cp.a = flux.a; %muon neutron pre-factor
cp.b = flux.b; %muon attenuation length
cp.negflux=flux.R; %flux negative muons
cp.totalflux=flux.phi; %total flux muons
%Also store the muons production rates from the code
cp.muon36(1,:)=cp.depthvector;

%store muon scaling factor
cp.SFmufast=flux.SFmufast;
cp.SFmuslow=flux.SFmuslow;

% Need to calculate the chemical compound factor
% Start by calculating the capture probability for target elements
% (numerator)

PcapturetargetK=cp.tni(1)* pp.table(9,12);
PcapturetargetCa=cp.tni(2)* pp.table(10,12);

%calculate the overall probability of capture (denominator)
Pcapturebulk=cp.ni.*pp.table(:,12);
%divide these two in the correct direction
cp.FcompoundK=PcapturetargetK/sum(Pcapturebulk);
cp.FcompoundCa=PcapturetargetCa/sum(Pcapturebulk);

% Find surface fluxes
cp.negflux0=interpolate(cp.depthvector,cp.negflux,0);
cp.totalflux0=interpolate(cp.depthvector,cp.totalflux,0);

% negative muon capture
P_negK = cp.negflux0.*pp.k_negpartial36K.*cp.FcompoundK.*pp.fstar36K;
P_negCa = cp.negflux0.*pp.k_negpartial36Ca.*cp.FcompoundCa.*pp.fstar36Ca;

%calculating fast muon contribution to chlorine (individually for Ca and K)
  z=cp.depthvector;

%see Marrero et al. 2016 on CRONUScalc for explanation of aalpha and Beta
aalpha = 1.0;
Beta = 1.0;

Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*0)) + 50.7.*(1-exp(-5.05e-7.*0)); %note, changed z to 0

P_fastK = cp.totalflux0.*Beta.*(Ebar.^aalpha).*pp.sigma036K.*cp.tni(1);
P_fastCa = cp.totalflux0.*Beta.*(Ebar.^aalpha).*pp.sigma036Ca.*cp.tni(2);
cp.P_fasttotal=P_fastK+P_fastCa;

%sum parts of the muon production

%Prodmu(jkl)=ProdmuK+ProdmuCa+P_negTi+P_negFe;
cp.muon36(1,1)=P_negCa+P_fastCa;
cp.muon36(2,1)=P_negK+P_fastK;

%
% NCa and NK are useful temporary variables.
%
NCa=cp.tni(2);
NK=cp.tni(1);
CU=sp36.ci(18);
CTh=sp36.ci(19);
NCl=cp.tni(5);
%
% Pnsf
%
cp.Pnsf=0.429*CU;
%
% X
%
% Depends on fi (cp.fi), Si (table(:,6), and Yn_iu (table(:,7)
%
cp.X=sum(pp.table(:,6).*cp.fi.*pp.table(:,7))/sum(pp.table(:,6).*cp.fi);
%
% Y
%
cp.Y=sum(pp.table(:,6).*cp.fi.*pp.table(:,8))/sum(pp.table(:,6).*cp.fi);
%
% Pnan
%
cp.Pnan=cp.X*CU+cp.Y*CTh;
%
% Pethr
%
cp.Pethr=(cp.Pnan+cp.Pnsf)*(1-cp.pEthss);
%
% Pthr
%
cp.Pthr=(cp.Pnan+cp.Pnsf)*cp.pEthss;
%
% N36r         Radiogenic production in atoms/gram.
%
cp.N36r=cp.Pthr*cp.fth/pp.lambda36Cl+cp.Pethr*cp.feth/pp.lambda36Cl;
%
% N36m         Measured Chlorine 36.
%
cp.N36m=sp36.concentration36;


