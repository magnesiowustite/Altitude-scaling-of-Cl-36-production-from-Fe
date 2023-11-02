%Code to calibrate the production rate of Cl-36 from K against Be-10 in qtz.
%for the Owens Valley samples and to propagate error.

function output = Owens_K()

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampledata = xlsread('Input Data.xlsx','Inputs','B2:M5');
sampledataK = xlsread('Input Data.xlsx','Inputs','B8:N11');

for i=1:3

data = sampledata(i,:);
dataK = sampledataK(i,:);

N10 = data(1,1);
eN10 = data(1,2);
thickness = data(1,7);
lambdaeff = data(1,8);
lambdaFe = data(1,9);
P10mu = data(1,10);
P10s = data(1,11);

N36K = dataK(1,1);
eN36K = dataK(1,2);
Misc_prod = dataK(1,12);
eMisc_prod = dataK(1,13);
Conc_K = dataK(1,5);
eConc_K = dataK(1,6);
S_K = dataK(1,7);
lambdaK = dataK(1,8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first compute the production rate
out = constant_exp_K(N10,N36K,thickness,Misc_prod,lambdaeff,lambdaK,P10mu,P10s,S_K,Conc_K);
out.P36m_K = out.P36m_Be;

%now compute partial derivatives and errors

%N10
dout = constant_exp_K((N10+N10*.001),N36K,thickness,Misc_prod,lambdaeff,lambdaK,P10mu,P10s,S_K,Conc_K);
dN10 = (dout.P36m_Be - out.P36m_K)/(N10*.001);
%N36m
dout = constant_exp_K(N10,(N36K+N36K*.001),thickness,Misc_prod,lambdaeff,lambdaK,P10mu,P10s,S_K,Conc_K);
dN36K = (dout.P36m_Be - out.P36m_K)/(N36K*.001);
%N36m_misc
dout = constant_exp_K(N10,N36K,thickness,(Misc_prod+Misc_prod*.001),lambdaeff,lambdaK,P10mu,P10s,S_K,Conc_K);
dN36misc = (dout.P36m_Be - out.P36m_K)/(Misc_prod*.001);
%Conc_K
dout = constant_exp_K(N10,N36K,thickness,Misc_prod,lambdaeff,lambdaK,P10mu,P10s,S_K,(Conc_K+Conc_K*.001));
dconcK = (dout.P36m_Be - out.P36m_K)/(Conc_K*.001);

%compute error
error10 = sqrt(dN10^2*eN10^2+dN36K^2*eN36K^2+dN36misc^2*eMisc_prod^2+dconcK^2*eConc_K^2);
output(i,:) = [out.P36m_K,error10];

end
