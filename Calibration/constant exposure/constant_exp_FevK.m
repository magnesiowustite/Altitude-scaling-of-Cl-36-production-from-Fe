%% Calibrates the Cl-36 production rate on K against Cl-36 in feldspar for a
%% constant-exposure geomorphic model. Does not consider differential attenuation.

function out = constant_exp_FevK(N36mt,N36fs,Prod_fs,Misc_mt,S_Fe,Conc_Fe,...
  thickness)

  lambdaeff = 176

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P36Fe = (N36mt/N36fs*(Prod_fs)-Misc_mt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate reference SLHL production ratios
out.P36m_Fe = P36Fe/S_Fe/Conc_Fe;
