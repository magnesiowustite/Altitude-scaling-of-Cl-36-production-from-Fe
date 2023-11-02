%Code to calibrate the production rate of Cl-36 from Fe against Be-10 in qtz.
%for the Owens Valley samples and to propagate error.

function output = Owens_FevBe()

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampledata = xlsread('Input Data.xlsx','Inputs','B2:O5');

for i=1:3

data = sampledata(i,:);

N10 = data(1,1);
eN10 = data(1,2);
N36m = data(1,3);
eN36m = data(1,4);
Misc_prod = data(1,5);
eMisc_prod = data(1,6);
thickness = data(1,7);
lambdaeff = data(1,8);
lambdaFe = data(1,9);
P10mu = data(1,10);
P10s = data(1,11);
S_Fe = data(1,12);
Conc_Fe = data(1,13);
eConc_Fe = data(1,14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first compute the production rate
out = constant_exp_FevBe(N10,N36m,thickness,Misc_prod,lambdaeff,lambdaFe,P10mu,P10s,S_Fe,Conc_Fe);
out.P36m_Be = out.P36m_Be;

%now compute partial derivatives and errors

%N10
dout = constant_exp_FevBe((N10+N10*.001),N36m,thickness,Misc_prod,lambdaeff,lambdaFe,P10mu,P10s,S_Fe,Conc_Fe);
dN10 = (dout.P36m_Be - out.P36m_Be)/(N10*.001);
%N36m
dout = constant_exp_FevBe(N10,(N36m+N36m*.001),thickness,Misc_prod,lambdaeff,lambdaFe,P10mu,P10s,S_Fe,Conc_Fe);
dN36m = (dout.P36m_Be - out.P36m_Be)/(N36m*.001);
%N36m_misc
dout = constant_exp_FevBe(N10,N36m,thickness,(Misc_prod+Misc_prod*.001),lambdaeff,lambdaFe,P10mu,P10s,S_Fe,Conc_Fe);
dN36misc = (dout.P36m_Be - out.P36m_Be)/(Misc_prod*.001);
%Conc_Fe
dout = constant_exp_FevBe(N10,N36m,thickness,Misc_prod,lambdaeff,lambdaFe,P10mu,P10s,S_Fe,(Conc_Fe+Conc_Fe*.001));
dconcFe = (dout.P36m_Be - out.P36m_Be)/(Conc_Fe*.001);

error10 = sqrt(dN10^2*eN10^2+dN36m^2*eN36m^2+dN36misc^2*eMisc_prod^2+dconcFe^2*eConc_Fe^2);
output(i,:) = [out.P36m_Be,error10];

end
