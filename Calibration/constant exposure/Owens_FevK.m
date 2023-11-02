%Code to calibrate the production rate of Cl-36 from Fe against Cl-36 from K
%for the Owens Valley samples and to propagate error.

function output = Owens_FevK()

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampledata = xlsread('Input Data.xlsx','Inputs','B2:P5');
sampledataK = xlsread('Input Data.xlsx','Inputs','B8:L11');

for i=1:3

data = sampledata(i,:);
dataK = sampledataK(i,:);

N36mt = data(1,3);
eN36mt = data(1,4);
Misc_mt = data(1,5);
eMisc_mt = data(1,6);
thickness = data(1,7);
lambdaeff = data(1,8);
lambdaFe = data(1,9);
S_Fe = data(1,12);
Conc_Fe = data(1,13);
eConc_Fe = data(1,14);

N36fs = dataK(1,1);
eN36fs = dataK(1,2);
Prod_fs = dataK(1,3);
eProd_fs = dataK(1,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first compute the production rate
out = constant_exp_FevK(N36mt,N36fs,Prod_fs,Misc_mt,S_Fe,Conc_Fe,...
   thickness,lambdaeff,lambdaFe);

%now compute partial derivatives and errors

%N36fs
dout = constant_exp_FevK(N36mt,(N36fs+0.001*N36fs),Prod_fs,Misc_mt,S_Fe,Conc_Fe,...
  thickness,lambdaeff,lambdaFe);
dN36fs = (dout.P36m_Fe - out.P36m_Fe)/(N36fs*.001);
%N36m
dout = constant_exp_FevK((N36mt+0.001*N36mt),N36fs,Prod_fs,Misc_mt,S_Fe,Conc_Fe,...
  thickness,lambdaeff,lambdaFe);
dN36mt = (dout.P36m_Fe - out.P36m_Fe)/(N36mt*.001);
%N36m_misc
dout = constant_exp_FevK(N36mt,N36fs,Prod_fs,(Misc_mt+0.001*Misc_mt),S_Fe,Conc_Fe,...
  thickness,lambdaeff,lambdaFe);
dmisc = (dout.P36m_Fe - out.P36m_Fe)/(Misc_mt*.001);
%Prod_fs
dout = constant_exp_FevK(N36mt,N36fs,(Prod_fs+0.001*Prod_fs),Misc_mt,S_Fe,Conc_Fe,...
  thickness,lambdaeff,lambdaFe);
dN36Prod_fs = (dout.P36m_Fe - out.P36m_Fe)/(Prod_fs*.001);
%Conc_Fe
dout = constant_exp_FevK(N36mt,N36fs,Prod_fs,Misc_mt,S_Fe,(Conc_Fe+0.001*Conc_Fe),...
  thickness,lambdaeff,lambdaFe);
dconcFe = (dout.P36m_Fe - out.P36m_Fe)/(Conc_Fe*.001);

%compute error
error = sqrt(dN36fs^2*eN36fs^2+dN36mt^2*eN36mt^2+dmisc^2*eMisc_mt^2+dconcFe^2*eConc_Fe^2);
output(i,:) = [out.P36m_Fe,error];

end
