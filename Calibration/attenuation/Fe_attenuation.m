% Calculates the Fe/Be attenuation length ratio in the atmosphere from 590 hPa
% to sea level at a cutoff rigidity of 6.3 GV.

function out = Fe_attenuation()

close all
load consts_LSD.mat;
w = 0.066;

for i = 590:10:1013.25

  out = LSDscaling_mod(i,6.3,462,w,consts,10);
  Be(1,44-(i-580)/10) = out;
  out = LSDscaling_mod(i,6.3,462,w,consts,56);
  Cl(1,44-(i-580)/10) = out;

end

  ratio = (Cl./Be);
  z = [1013.25:-10:590];

##  hold all
##  plot(ratio,z)
##  set (gca (), "ydir", "reverse")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data and model for parameterization for Fe
X = z'; %setup fitting depths
P_Fe = Cl';

leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2)); %function to fit
F = leasqrfunc;

%Generate Fit
pin = [0,160]; % initial value estimates
[f,p,r21]  = leasqr(X,P_Fe,pin,F); %least-squares fit

%Outputs
a_Fe = p(1,1); %pre-factor
b_Fe = p(2,1); %attention length
r21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data and model for parameterization for Be
X = z'; %setup fitting depths
P_Be = Be';

leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2)); %function to fit
F = leasqrfunc;

%Generate Fit
pin = [0,160]; % initial value estimates
[f,p,r21]  = leasqr(X,P_Be,pin,F); %least-squares fit

%Outputs
a_Be = p(1,1); %pre-factor
b_Be = p(2,1); %attention length
r21;

%Calculate attenuation length ratio
out = b_Fe./b_Be;

%Plot best-fit model
P_Fe = a_Fe*e.^(-z./b_Fe);
P_Be = a_Be*e.^(-z./b_Be);













