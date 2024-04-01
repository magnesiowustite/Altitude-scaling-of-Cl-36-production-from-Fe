%% Calibrates the Cl-36 production rate from low-energy neutrons against Be-10...
%% in qtz for a constant-exposure geomorphic model. Does not consider
%% differential attenuation

function out = constant_exp_K(N10,N36m,thickness,Misc_prod,lambdaeff,lambdaK,...
  Pmu,P10s,S_K,Conc_K)

%Decay constants
lambda10 = -log(1/2)/1.39E12;
lambda36 = -log(1/2)/3.01E12;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute thickness factors
F_Be = lambdaeff/(thickness)*(1-e^(-thickness/lambdaeff));
F_K = lambdaeff/(thickness)*(1-e^(-thickness/lambdaeff));

%Solve Be-10 for exposure age
t = 10;
N10n = 0;
while (abs(N10-N10n) >= 1);
N10n = F_Be*P10s/lambda10*(1-e^(-lambda10*t))+Pmu/lambda10*(1-e^(-lambda10*t));
t = N10/N10n*t; %note, we assume muon attenuation is neglible across the sample thickness
endwhile

P36Be = N36m/F_K*lambda36/(1-e^(-lambda36*t)); %Solve for Cl-36 production rate in mt.
P36Be = P36Be - Misc_prod - S_K*Conc_K*155.7; %subtract misc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out.P36m_Be = P36Be;
