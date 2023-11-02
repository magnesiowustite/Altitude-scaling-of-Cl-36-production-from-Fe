% Code for calculating muon production parameters for Cl-36
% modified from Marrero et al. (2016)
% - Muonfluxsato changed to 3-exponential term muon model and
% - flux depth vector changed to flux0 for muon production on Ca and K.

function cp=comppars36_muons(pp,sp36)

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
% Do I need to add something here?

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

flux=Muon_model_muons(cp.depthvector,sp36.P,RcEst,pp.SPhiInf);
%store the output fluxes that we need
%save coefficients and attenuation lengths
cp.a = flux.a;
cp.b = flux.b;
cp.c = flux.c;
cp.d = flux.d;
cp.e = flux.e;
cp.f = flux.f;

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

% negative muon capture surface production rates
cp.P_negK_a = cp.a.*pp.k_negpartial36K.*cp.FcompoundK.*pp.fstar36K;
cp.P_negK_c = cp.c.*pp.k_negpartial36K.*cp.FcompoundK.*pp.fstar36K;
cp.P_negCa_a = cp.a.*pp.k_negpartial36Ca.*cp.FcompoundCa.*pp.fstar36Ca;
cp.P_negCa_c = cp.c.*pp.k_negpartial36Ca.*cp.FcompoundCa.*pp.fstar36Ca;
cp.P_negtotal36_a = cp.P_negK_a+cp.P_negCa_a;
cp.P_negtotal36_c = cp.P_negK_c+cp.P_negCa_c;

%calculating fast muon contribution to chlorine (individually for Ca and K)
  z=cp.depthvector;

%see Marrero et al. 2016 on CRONUScalc for explanation of aalpha and Beta
aalpha = 1.0;
Beta = 1.0;

Ebar = 7.6 + 321.7.*(1 - exp(-8.059e-6.*0)) + 50.7.*(1-exp(-5.05e-7.*0)); %note, changed z to 0

%Fast muon production rates
cp.P_fastK = cp.e.*Beta.*(Ebar.^aalpha).*pp.sigma036K.*cp.tni(1);
cp.P_fastCa = cp.e.*Beta.*(Ebar.^aalpha).*pp.sigma036Ca.*cp.tni(2);
cp.P_fasttotal36 = cp.P_fastK+cp.P_fastCa;

