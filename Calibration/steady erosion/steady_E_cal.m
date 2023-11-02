%% This code calibrates the production rate of Cl-36 from Fe against Be-10 in
%% qtz and Cl-36 in feldspar and the Cl-36 production rate from K against
%% Be-10 in quartz at Mt. Evans using a steady-state erosion geomophic model.
%% Includes reaction-specific spallation attenuation lengths.

function out = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,...
  misc_fs_sp,misc_mt_sp,Conc_Fe)

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%muons
muons = muon_parameters();
P36mus_a = muons.P_neg36_a;
P36mus_c = muons.P_neg36_c;
P10mus_a = muons.P_neg10_a;
P10mus_c = muons.P_neg10_c;
Lmu_b = muons.b;
Lmu_d = muons.d;
P10muf = muons.P_fast10;
P36muf = muons.P_fast36;
Lmu_f = muons.f;

%neutrons
sampledata = xlsread('Mt Evans.xlsx','Inputs','B2:AM49');
sf = xlsread('Mt Evans.xlsx','Scaling','B2:F4');
pp = physpars_NC();
sp36 = samppars36_NC(sampledata);
cp36 = comppars36_NC(pp,sp36);
out = prodz36_all(pp,sp36,cp36,sf);

Lfe = out.Lambdafe;
Lmu = out.Lambdamu;
Lth = out.Lth;
Leth = out.Leth;
k1 = out.k1;
k2 = out.k2;
k3 = out.k3;
k4 = out.k4;
k5 = out.k5;
k6 = out.k6;
k7 = out.k7;

%Correct feldspar Cl-36 conc. for radiogenic neutron production
N36F = N36F-out.R;

%reaction-specific attenuation lengths
[ratio] = Argento_attenuation;

Lef = Lfe;
lambdaK = ratio*Lef;
lambdaFe = Lef*0.91;
lambda10 = -log(1/2)/1.389E6;
lambda36 = -log(1/2)/3.01E5;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute the thickness factors for Be
F = Lef/(thickness)*(1-e^(-thickness/Lef));
F_Fe = lambdaFe/(thickness)*(1-e^(-thickness/lambdaFe));
Fmu_a = Lmu_b/(thickness)*(1-e^(-thickness/Lmu_b));
Fmu_c = Lmu_d/(thickness)*(1-e^(-thickness/Lmu_d));
Fmu_f = Lmu_f/(thickness)*(1-e^(-thickness/Lmu_f));

%Solve Be-10 for erosion rate
E = .0001;
N10n = 0;
while (abs(N10-N10n) >= 1)
N10n = F*P10s/(lambda10+E/Lef) + Fmu_a*P10mus_a/(lambda10+E/Lmu_b)...
+Fmu_c*P10mus_c/(lambda10+E/Lmu_d)+Fmu_f*P10muf/(lambda10+E/Lmu_f);
E = N10n/N10*E;
endwhile

%Solve for Cl-36 production rate in mt. and subtact misc spall. in magnetite
P36Be = (N36m-F*misc_mt_sp/(lambda36+E/Lef))*(lambda36+E/lambdaFe)/F_Fe;

%K production rate calibration
F_K = lambdaK/(thickness)*(1-e^(-thickness/lambdaK));
F_eth = Leth/(thickness)*(1-e^(-thickness/Leth));
F_th = Lth/(thickness)*(1-e^(-thickness/Lth));
F_mu = Lmu/(thickness)*(1-e^(-thickness/Lmu));

%Solve for K production rate
P36KBe = (N36F-(F*misc_fs_sp/(lambda36+E/Lef)+Fmu_a*P36mus_a/(lambda36+E/Lmu_b)...
+Fmu_c*P36mus_c/(lambda36+E/Lmu_d)+Fmu_f*P36muf/(lambda36+E/(Lmu_f))+F*k1/(lambda36+E/(Lef))...
+F_eth*k2/(lambda36+E/(Leth))+F_mu*k3/(lambda36+E/(Lmu))+F*k4/(lambda36+E/(Lef))+...
F_eth*k5/(lambda36+E/(Leth))+F_th*k6/(lambda36+E/(Lth))+F_mu*k7/(lambda36+E/(Lmu))))/F_K*...
(lambda36+E/lambdaK);

%Solve Cl-36 for erosion rate
E = .0001;
N36Fn = 0;
while (abs(N36F-N36Fn) >= 1)
N36Fn = F_K*P36s/(lambda36+E/lambdaK)+ F*misc_fs_sp/(lambda36+E/Lef) + Fmu_a*P36mus_a/(lambda36+E/Lmu_b)...
+Fmu_c*P36mus_c/(lambda36+E/Lmu_d)+Fmu_f*P36muf/(lambda36+E/(Lmu_f))+F*k1/(lambda36+E/(Lef))...
+F_eth*k2/(lambda36+E/(Leth))+F_mu*k3/(lambda36+E/(Lmu))+F*k4/(lambda36+E/(Lef))...
+F_eth*k5/(lambda36+E/(Leth))+F_th*k6/(lambda36+E/(Lth))+F_mu*k7/(lambda36+E/(Lmu));
E = N36Fn/N36F*E;
endwhile

%Solve for Cl-36 production rate in mt.
%subtract misc spall. in magnetite
P36K = (N36m-F*misc_mt_sp/(lambda36+E/Lef))*(lambda36+E/lambdaFe)/F_Fe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate reference SLHL production ratios
out.P36m_Be = P36Be/Conc_Fe/S_Fe;
out.P36m_K = P36K/Conc_Fe/S_Fe;
out.P36K_Be = P36KBe/Conc_K/S_K;
