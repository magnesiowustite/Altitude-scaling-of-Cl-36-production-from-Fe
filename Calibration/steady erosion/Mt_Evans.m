%% This code calculates production rates from the Mt. Evans sample and
%% propagates errors by a finite difference approximation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampledata = xlsread('Input data.xlsx','Inputs','B5:P5');
sampledataK = xlsread('Input data.xlsx','Inputs','B11:M11');

N10 = sampledata(1,1); %conc. Be-10
eN10 = sampledata(1,2); %error conc. Be-10
N36m = sampledata(1,3); %conc. Cl-36 in mt.
eN36m = sampledata(1,4); %error Cl-36 in mt.
thickness = sampledata(1,7); %thickness
P10s = sampledata(1,11); %spallation production of Be-10
S_Fe = sampledata(1,12); %spallation scaling factor for Fe
Conc_Fe = sampledata(1,13); %concentration of Fe
eConc_Fe = sampledata(1,14); %error on concentration of Fe
misc_mt_sp = sampledata(1,15); %misc. spallation production in magnetite

N36F = sampledataK(1,1); %conc. Cl-36 in fs.
eN36F = sampledataK(1,2); %error Cl-36 in fs.
misc_fs_sp = sampledataK(1,11); %misc. spallation production in feldspar
P36s = sampledataK(1,9); %spallation production in feldspar
eP36s = sampledataK(1,10); %error spallation production in feldspar
Conc_K = sampledataK(1,5); %concentration of K
eConc_K = sampledataK(1,6); %error on concentration of K
S_K = sampledataK(1,7); %spallation scaling factor for K

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%first compute production rates
out = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,misc_fs_sp,...
misc_mt_sp,Conc_Fe);

%now compute partial derivatives and errors
%Fe, Be-10 normalized
%N10
dout = steady_E_cal((N10+N10*.001),P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN10 = (dout.P36m_Be - out.P36m_Be)/(N10*.001);
%N36m
dout = steady_E_cal(N10,P10s,S_Fe,(N36m+N36m*.001),N36F,Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN36m = (dout.P36m_Be - out.P36m_Be)/(N36m*.001);
%Conc_Fe
dout = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,misc_fs_sp,...
misc_mt_sp,(Conc_Fe+Conc_Fe*.001));
dconcFe = (dout.P36m_Be - out.P36m_Be)/(Conc_Fe*.001);

error10 = sqrt(dN10^2*eN10^2+dN36m^2*eN36m^2+dconcFe^2*eConc_Fe^2);

%Fe, Cl-36 normalized
%N36m
dout = steady_E_cal(N10,P10s,S_Fe,(N36m+N36m*.001),N36F,Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN36m = (dout.P36m_K - out.P36m_K)/(N36m*.001);
%N36F
dout = steady_E_cal(N10,P10s,S_Fe,N36m,(N36F+N36F*.001),Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN36F = (dout.P36m_K - out.P36m_K)/(N36F*.001);
%N36F_target
dout = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,Conc_K,S_K,(P36s+P36s*.001),thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN36s = (dout.P36m_K - out.P36m_K)/(P36s*.001);
%Conc_Fe
dout = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,misc_fs_sp,...
misc_mt_sp,(Conc_Fe+Conc_Fe*.001));
dconcFe = (dout.P36m_K - out.P36m_K)/(Conc_Fe*.001);

error36Fe = sqrt(dN36F^2*eN36F^2+dN36m^2*eN36m^2+dN36s^2*eP36s^2+dconcFe^2*eConc_Fe^2);

%K, Be-10 normalized
%N10
dout = steady_E_cal((N10+N10*.001),P10s,S_Fe,N36m,N36F,Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN10 = (dout.P36K_Be - out.P36K_Be)/(N10*.001);
%N36m
dout = steady_E_cal(N10,P10s,S_Fe,N36m,(N36F+N36F*.001),Conc_K,S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dN36m = (dout.P36K_Be - out.P36K_Be)/(N36F*.001);
%Conc_K
dout = steady_E_cal(N10,P10s,S_Fe,N36m,N36F,(Conc_K+Conc_K*.001),S_K,P36s,thickness,...
misc_fs_sp,misc_mt_sp,Conc_Fe);
dconcK = (dout.P36K_Be - out.P36K_Be)/(Conc_K*.001);

error36K = sqrt(dN10^2*eN10^2+dN36m^2*eN36m^2+dconcK^2*eConc_K^2);

output = [out.P36m_K,error36Fe,out.P36m_Be,error10,out.P36K_Be,error36K]
