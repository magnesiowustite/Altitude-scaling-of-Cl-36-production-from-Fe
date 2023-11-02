%% Calibrates the Cl-36 production rate on Fe against Be-10 in qtz for a
%% constant-exposure geomorphic model.

function out = constant_exp_FevBe(N10,N36m,thickness,Misc_prod,lambdaeff,...
  lambdaFe,Pmu,P10s,S_Fe,Conc_Fe)

%Decay constants
lambda10 = -log(1/2)/1.39E6;
lambda36 = -log(1/2)/3.01E5;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute thickness factors
F_Be = lambdaeff/(thickness)*(1-e^(-thickness/lambdaeff));
F_Fe = lambdaFe/(thickness)*(1-e^(-thickness/lambdaFe));

%Solve Be-10 for erosion rate
t = 10;
N10n = 0;
while (abs(N10-N10n) >= 1);
N10n = F_Be*P10s/lambda10*(1-e^(-lambda10*t))+Pmu/lambda10*(1-e^(-lambda10*t));
t = N10/N10n*t; %note, we assume muon attenuation is neglible across the sample thickness
endwhile

P36Be = N36m/F_Fe*lambda36/(1-e^(-lambda36*t)); %Solve for the Cl-36 production rate in mt.
P36Be = P36Be - Misc_prod; %subtract misc. production pathways

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate reference SLHL production ratios by normalizing by the scaling factor
%and Fe concentration
out.P36m_Be = P36Be/S_Fe/Conc_Fe;
