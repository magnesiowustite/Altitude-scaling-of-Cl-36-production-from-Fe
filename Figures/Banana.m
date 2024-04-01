%% Plot Banana at Mt. Evans %%

%Makes figure 5, a 2-isotope Cl-36/Be-10 plot for the Mt. Evans sample.

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setup Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sampledata = xlsread('Input Data.xlsx','Inputs','B5:P5');
sampledataK = xlsread('Input Data.xlsx','Inputs','B11:M11');

N10 = sampledata(1,1); %conc. Be-10
eN10 = sampledata(1,2); %error conc. Be-10
N36m = sampledata(1,3); %conc. Cl-36 in mt.
eN36m = sampledata(1,4); %error Cl-36 in mt.
N36F = sampledataK(1,1); %conc. Cl-36 in fs.
eN36F = sampledataK(1,2); %error Cl-36 in fs.
thickness = sampledata(1,7); %thickness
misc_mt_sp = sampledata(1,15); %misc. spallation production in magnetite
misc_fs_sp = sampledataK(1,11); %misc. spallation production in feldspar
P36s = sampledataK(1,9); %spallation production in feldspar
eP36s = sampledataK(1,10); %error spallation production in feldspar
P10s = sampledata(1,11); %spallation production of Be-10
S_Fe = sampledata(1,12); %spallation scaling factor for Fe
Conc_K = sampledataK(1,5); %concentration of K
eConc_K = sampledataK(1,6); %error on concentration of K
S_K = sampledataK(1,7); %spallation scaling factor for K
Conc_Fe = sampledata(1,13); %concentration of Fe
eConc_Fe = sampledata(1,14); %error on concentration of Fe

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
sampledata = xlsread('Mt Evans.xlsx','Inputs','B2:AM2');
sf         = xlsread('Mt Evans.xlsx','Scaling','B2:F2');
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

%calculate attenuation length ratios
[ratio_K] = Argento_attenuation;

Lef = Lfe; %Be-10 attenuation length
lambdaK = ratio_K*Lef;

%Decay constants
lambda10 = -log(1/2)/1.389E6;
lambda36 = -log(1/2)/3.01E5;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Steady State Erosion

%Compute thickness factor for Be
F = Lef/(thickness)*(1-e^(-thickness/Lef));
Fmu_a = Lmu_b/(thickness)*(1-e^(-thickness/Lmu_b));
Fmu_c = Lmu_d/(thickness)*(1-e^(-thickness/Lmu_d));
Fmu_f = Lmu_f/(thickness)*(1-e^(-thickness/Lmu_f));

%Be-10 for erosion rate
for i = 1000:1000:100E7;
E = 1E-8*(i);
N10n(i/1000) = F*P10s/(lambda10+E/Lef) + Fmu_a*P10mus_a/(lambda10+E/Lmu_b)...
+Fmu_c*P10mus_c/(lambda10+E/Lmu_d)+Fmu_f*P10muf/(lambda10+E/Lmu_f);
end

%K production rate calibration
F_K = lambdaK/(thickness)*(1-e^(-thickness/lambdaK));
F_eth = Leth/(thickness)*(1-e^(-thickness/Leth));
F_th = Lth/(thickness)*(1-e^(-thickness/Lth));
F_mu = Lmu/(thickness)*(1-e^(-thickness/Lmu));

%Cl-36 for erosion rate
for i = 1000:1000:100E7;
E = 1E-8*(i);
N36Fn(i/1000) = F_K*P36s/(lambda36+E/lambdaK)+ F*misc_fs_sp/(lambda36+E/Lef) + Fmu_a*P36mus_a/(lambda36+E/Lmu_b)...
+Fmu_c*P36mus_c/(lambda36+E/Lmu_d)+Fmu_f*P36muf/(lambda36+E/(Lmu_f))+F*k1/(lambda36+E/(Lef))...
+F_eth*k2/(lambda36+E/(Leth))+F_mu*k3/(lambda36+E/(Lmu))+F*k4/(lambda36+E/(Lef))...
+F_eth*k5/(lambda36+E/(Leth))+F_th*k6/(lambda36+E/(Lth))+F_mu*k7/(lambda36+E/(Lmu));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

ratio_E = N36Fn./N10n;

%Constant exposure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute thickness factor
F_Be = F;

%Solve Be-10 for exposure age
for i = 1000:1000:10E6;
t = i;
N10nx(i/1000) = F*P10s/lambda10*(1-e^(-lambda10*t))...
+ Fmu_a*P10mus_a/lambda10*(1-e^(-lambda10*t))...
+Fmu_c*P10mus_c/lambda10*(1-e^(-lambda10*t))...
+Fmu_f*P10muf/lambda10*(1-e^(-lambda10*t));
end

%Solve Cl-36 for exposure age
for i = 1000:1000:10E6;
t = i;
N36Fnx(i/1000) = F_K*P36s/lambda36*(1-e^(-lambda36*t))+...
F*misc_fs_sp/lambda36*(1-e^(-lambda36*t)) + Fmu_a*P36mus_a/lambda36*(1-e^(-lambda36*t))...
+Fmu_c*P36mus_c/lambda36*(1-e^(-lambda36*t))+Fmu_f*P36muf/lambda36*(1-e^(-lambda36*t))...
+F*k1/lambda36*(1-e^(-lambda36*t))+F_eth*k2/lambda36*(1-e^(-lambda36*t))...
+F_mu*k3/lambda36*(1-e^(-lambda36*t))+F*k4/lambda36*(1-e^(-lambda36*t))...
+F_eth*k5/lambda36*(1-e^(-lambda36*t))+F_th*k6/lambda36*(1-e^(-lambda36*t))...
+F_mu*k7/lambda36*(1-e^(-lambda36*t));
end

ratio_exp = N36Fnx./N10nx;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
h1 = semilogx(N10n,ratio_E,'color','k','linewidth',1.5,'linestyle','-','marker','none');
h2 = semilogx(N10nx,ratio_exp,'color','k','linestyle','--','linewidth',1.5,'marker','none');

set(gca,"fontsize", 16,'fontweight','normal', "linewidth", 2, 'XMinorTick', 'on',...
'YMinorTick','on','TickLength',[0.025, 0.0025]);

box on
xlim([1e5,1e9/2]);
ylim([0,5]);

%data
N10 = xlsread('Input data.xlsx','Inputs','B5');
N10e = xlsread('Input data.xlsx','Inputs','C5');
N36 = xlsread('Input data.xlsx','Inputs','B11');
N36e = xlsread('Input data.xlsx','Inputs','C11');
ratio = N36/N10;
ratioe = ratio*sqrt((N10e/N10)^2+(N36e/N36)^2);

h3 = errorbar(N10, ratio, 2*N10e, 2*ratioe,"#~>-*");
set(h3,'color','k','linestyle','none','linewidth',1.5,'marker','o',...
'markerfacecolor','r','markersize',10,...
'markeredgecolor','k');

ylabel('^{36}Cl_{K-fs.}/^{10}Be_{qtz}',"fontweight","normal","fontsize",16);
xlabel('^{10}Be_{qtz} (atoms g^{-1})',"fontweight","normal","fontsize",16);

text(150000,4.55,'steady erosion','fontsize',16)
text(150000,3.7,'constant exposure','fontsize',16)
text(2500000,4.3,'Mt. Evans','fontsize',16)

  % Write results to file
    set(gcf, "paperunits", "points", "papersize", [900, 600], 'PaperPosition', [0 0 900 600]);
   print -dpdf -color Banana.pdf
   print -dpng -color Banana.png


