function [] = Fluxes_Plot()

%Makes Figure 1, a plot of normalized particles fluxes and excitation functions
clf

% Some parameters, note that Rc and SPhi are averages for last 21000 years
% at the Sierra Nevada locations.

  Rc = 6.3 % cutoff rigidity
  s = 462 % SPhi
  w = 0.066 % water content of the ground
  nuclide = 56

load consts_LSD.mat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Neutrons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define pressures that we will plot
pressures = [1033,830.9,683.2,603.0]
%setup subplot for neutrons
subplot(2,1,1);

h = subplot(2,1,1);
p = get(h, 'position')
p(3) = .8
p(4) = .2*1200/900
set(h, 'position', p);

for i = 1:4

  h = pressures(1,i)

% The next section is from Lifton et al. (2014):

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sato et al. (2008) Neutron Spectrum
% Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, 2013
% Purdue University, nlifton@purdue.edu

% Copyright 2013, Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).

x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)

% E = logspace(-8,5,1000);
E = logspace(0,5.3010,200);
% E = [1.1295 11.295 112.95 1129.5 11295];

% Flatten low rigidities.

lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

%nflux = zeros(length(Rc));

%w = 0.2; % water content, from 0-1
%s = 1700;
%Rc = 12;
%x = 1030;
Et = 2.5e-8; % Thermal Neutron Energy in MeV

% Integrated neutron flux <15 MeV

smin = 400; %units of MV
smax = 1200; %units of MV

a6 = 1.8882e-4;
a7 = 4.4791e-1;
a8 = 1.4361e-3;
a12 = 1.4109e-2;

b11min = 2.5702e1;
b11max = -6.9221;
b12min = -5.0931e-1;
b12max = 1.1336;
b13min= 7.4650;
b13max = 2.6961e1;
b14min = 1.2313e1;
b14max = 1.1746e1;
b15min = 1.0498;
b15max = 2.6171;
b21min = 6.5143e-3;
b21max = 5.3211e-3;
b22min = 3.3511e-5;
b22max = 8.4899e-5;
b23min = 9.4415e-4;
b23max = 2.0704e-3;
b24min = 1.2088e1;
b24max = 1.1714e1;
b25min = 2.7782;
b25max = 3.8051;
b31min = 9.8364e-1;
b31max = 9.7536e-1;
b32min = 1.4964e-4;
b32max = 6.4733e-4;
b33min = -7.3249e-1;
b33max = -2.2750e-1;
b34min = -1.4381;
b34max = 2.0713;
b35min = 2.7448;
b35max = 2.1689;
b41min = 8.8681e-3;
b41max = 9.1435e-3;
b42min = -4.3322e-5;
b42max = -6.4855e-5;
b43min = 1.7293e-2;
b43max = 5.8179e-3;
b44min = -1.0836;
b44max = 1.0168;
b45min = 2.6602;
b45max = 2.4504;

b121 = 9.31e-1;
b122 = 3.70e-2;
b123 = -2.02;
b124 = 2.12;
b125 = 5.34;
b131 = 6.67e-4;
b132 = -1.19e-5;
b133 = 1.00e-4;
b134 = 1.45;
b135 = 4.29;

% Basic Spectrum

b51 = 9.7337e-4;
b52 = -9.6635e-5;
b53 = 1.2121e-2;
b54 = 7.1726;
b55 = 1.4601;
b91 = 5.7199e2;
b92 = 7.1293;
b93 = -1.0703e2;
b94 = 1.8538;
b95 = 1.2142;
b101 = 6.8552e-4;
b102 = 2.7136e-5;
b103 = 5.7823e-4;
b104 = 8.8534;
b105 = 3.6417;
b111 = -5.0800e-1;
b112 = 1.4783e-1;
b113 = 1.0068;
b114 = 9.1556;
b115 = 1.6369;

c1 = 2.3555e-1; % lethargy^-1
c2 = 2.3779; % MeV
c3 = 7.2597e-1;
c5 = 1.2391e2; % MeV
c6 = 2.2318; % MeV
c7 = 1.0791e-3; % lethargy^-1
c8 = 3.6435e-12; % MeV
c9 = 1.6595;
c10 = 8.4782e-8; % MeV
c11 = 1.5054;

% Ground-Level Spectrum

h31 = -2.5184e1;
h32 = 2.7298;
h33 = 7.1526e-2;
h51 = 3.4790e-1;
h52 = 3.3493;
h53 = -1.5744;

g1 = -0.023499;
g2 = -0.012938;
g3 = 10.^(h31 + h32./(w + h33));
g4 = 9.6889e-1;
g5 = h51 + h52.*w + h53.*(w.^2);

fG = 10.^(g1 + g2.*log10(E./g3).*(1-tanh(g4.*log10(E./g5))));

% Thermal Neutron Spectrum

h61 = 1.1800e-1;
h62 = 1.4438e-1;
h63 = 3.8733;
h64 = 6.5298e-1;
h65 = 4.2752e1;

g6 = (h61 + h62.*exp(-h63.*w))./(1 + h64.*exp(-h65.*w));

PhiT = g6.*((E./Et).^2).*exp(-E./Et);

% Total Ground-Level Flux

PhiB = zeros(1,length(E));
PhiG = zeros(1,length(E));
PhiGMev = zeros(1,length(E));
p10n = zeros(1,length(E));
p14n = zeros(1,length(E));
p26n = zeros(1,length(E));
p3n = zeros(1,length(E));
p36Can = zeros(1,length(E));
p36Kn = zeros(1,length(E));
p36Tin = zeros(1,length(E));
p36Fen = zeros(1,length(E));

for a = 1:length(Rc)

    a1min = b11min + b12min.*Rc(a) + b13min./(1 + exp((Rc(a) - b14min)./b15min));
    a1max = b11max + b12max.*Rc(a) + b13max./(1 + exp((Rc(a) - b14max)./b15max));
    a2min = b21min + b22min.*Rc(a) + b23min./(1 + exp((Rc(a) - b24min)./b25min));
    a2max = b21max + b22max.*Rc(a) + b23max./(1 + exp((Rc(a) - b24max)./b25max));
    a3min = b31min + b32min.*Rc(a) + b33min./(1 + exp((Rc(a) - b34min)./b35min));
    a3max = b31max + b32max.*Rc(a) + b33max./(1 + exp((Rc(a) - b34max)./b35max));
    a4min = b41min + b42min.*Rc(a) + b43min./(1 + exp((Rc(a) - b44min)./b45min));
    a4max = b41max + b42max.*Rc(a) + b43max./(1 + exp((Rc(a) - b44max)./b45max));

    a5 = b51 + b52.*Rc(a) + b53./(1 + exp((Rc(a) - b54)./b55));
    a9 = b91 + b92.*Rc(a) + b93./(1 + exp((Rc(a) - b94)./b95));
    a10 = b101 + b102.*Rc(a) + b103./(1 + exp((Rc(a) - b104)./b105));
    a11 = b111 + b112.*Rc(a) + b113./(1 + exp((Rc(a) - b114)./b115));

    b5 = b121 + b122.*Rc(a) + b123./(1 + exp((Rc(a) - b124)./b125));
    b6 = b131 + b132.*Rc(a) + b133./(1 + exp((Rc(a) - b134)./b135));

    c4 = a5 + a6.*x./(1 + a7.*exp(a8.*x)); % lethargy^-1
    c12 = a9.*(exp(-a10.*x) + a11.*exp(-a12.*x)); % MeV

    PhiLmin = a1min.*(exp(-a2min.*x) - a3min.*exp(-a4min.*x)); %length of Rc
    PhiLmax = a1max.*(exp(-a2max.*x) - a3max.*exp(-a4max.*x)); %length of Rc

    f3 = b5 + b6.*x;
    f2 = (PhiLmin - PhiLmax)./(smin.^f3 - smax.^f3);
    f1 = PhiLmin - f2.*smin.^f3;

    PhiL = f1 + f2.*s(a).^f3;

    PhiB = (c1.*(E./c2).^c3).*exp(-E./c2) + c4.*exp((-(log10(E) - log10(c5)).^2)./(2.*(log10(c6)).^2))...
    + c7.*log10(E./c8).*(1 + tanh(c9.*log10(E./c10))).*(1 - tanh(c11.*log10(E./c12)));

    PhiG = PhiL.*(PhiB.*fG + PhiT);
    PhiGMev = PhiG./E;

    clipindex = find(E <= 1, 1, 'last' ); %Make sure the clip index is consistent with the definition of E above

end

N.E = E;
N.PhiG(i,:) = PhiG./PhiG(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the neutron spectra and excitation functions
hold on
box on

  [AX1,H1,H2] = plotyy(E,PhiG./PhiG(1,1)./N.PhiG(1,:),E,consts.FenxCl36,@semilogx,@semilogy);
 set(AX1(2),'fontsize',18,'ycolor','k')
 set(H2,'linestyle',':','linewidth',2,'color','b','linestyle','--')
 set(H1,'linewidth',2,'color',[.5,.5,.5])

 [AX2,H1,H5] = plotyy(E,PhiG./PhiG(1,1)./N.PhiG(1,:),E,consts.O16nxBe10,@semilogx,@semilogx);
 set(H1,'linewidth',2,'color','k')
 set(H5,'linewidth',2,'color','k','linestyle','-.')

  [AX2,H3,H4] = plotyy(E,PhiG./PhiG(1,1)./N.PhiG(1,:),E,consts.KnxCl36,@semilogx,@semilogx);
 set(AX2(1),'ycolor',[0,0,0],'fontsize',18,'XMinorTick','on','linewidth',3)
 set(AX2(2),'ycolor','k','fontsize',18,'XMinorTick','on','linewidth',3)
 set(H4,'linestyle',':','linewidth',2,'color',[.9,.5,.2])
 set(H3,'linewidth',2,'color',i*[.2,.2,.2])
 set(AX2(1),'YLim',[1 2.2])
 set(AX2(2),'YLim',[.1 1000])
 set(AX2,'XLim',[0 100000])

 xlabel('Energy (MeV)','fontsize',20)
 ylabel(AX1(1),'Normalized Neutron Flux')
 ylabel(AX1(2),'millibarns')
 legend([H2,H5,H4,H1],'Fe(n,x)^{36}Cl','^{16}O(n,x)^{10}Be','K(n,x)^{36}Cl','normalized flux','Location','southwest')
 legend boxoff
 text(20000,2.08,'4300 m','fontsize',18)
 text(20000,1.85,'3300 m','fontsize',18)
 text(20000,1.43,'1700 m','fontsize',18)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Protons %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define pressures
pressures = [1033,830.9,683.2,603.0]
%Setup second subplot
subplot(2,1,2);

%Adjust position of subplots
h = subplot(2,1,2);
p = get(h, 'position')
p(2) = p(2) + 0.3;
p(3) = .8
p(4) = .05*1200/775
set(h, 'position', p);

%loop over pressures
for i = 1:4

  h = pressures(1,i)

% Sato et al. (2008) Neutron Spectrum
% Analytical Function Approximation (PARMA)
% Implemented in MATLAB by Nat Lifton, 2013
% Purdue University, nlifton@purdue.edu

% Copyright 2013, Purdue University
% All rights reserved
% Developed in part with funding from the National Science Foundation.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License, version 3,
% as published by the Free Software Foundation (www.fsf.org).

x = h.*1.019716; % Convert pressure (hPa) to atm depth (g/cm2)

E = logspace(0,5.3010,200);
% E = [1.1295 11.295 112.95 1129.5 11295];

A = 1;
Z = 1;
Ep = 938.27; % Rest mass of a proton
U = (4-1.675).*pi.*A./Z.*1e-7; % Unit conversion factor

% Flatten low rigidities.

lowRc = find(Rc < 1.0);
Rc(lowRc) = 1.0 + zeros(size(lowRc));

% Primary spectrum

smin = 400; %units of MV
smax = 1200; %units of MV

a1 = 2.1153;
a2 = 4.4511e-1;
a3 = 1.0064e-2;
a4 = 3.9564e-2;
a5 = 2.9236;
a6 = 2.7076;
a7 = 1.2663e4;
a8 = 4.8288e3;
a9 = 3.2822e4;
a10 = 7.4378e3;
a11 = 3.4643;
a12 = 1.6752;
a13 = 1.3691;
a14 = 2.0665;
a15 = 1.0833e2;
a16 = 2.3013e3;

Etoa = E + a1.*x;
Rtoa = 0.001.*sqrt((A.*Etoa).^2 + 2.*A.*Ep.*Etoa)./Z;

Elis = zeros(1,length(E));
Beta = zeros(1,length(E));
Rlis = zeros(1,length(E));
phiTOA = zeros(1,length(E));
phiLIS = zeros(1,length(E));
phiSec = zeros(1,length(E));
phiPtot = zeros(1,length(E));
p10p = zeros(1,length(E));

% Secondary Spectrum

c11 = 1.2560;
c12 = 3.2260e-3;
c13 = -4.8077e-6;
c14 = 2.2825e-9;
c21 = 4.3783e-1;
c22 = -5.5759e-4;
c23 = 7.8388e-7;
c24 = -3.8671e-10;
c31 = 1.8102e-4;
c32 = -5.1754e-7;
c33 = 7.5876e-10;
c34 = -3.8220e-13;
c41 = 1.7065;
c42 = 7.1608e-4;
c43 = -9.3220e-7;
c44 = 5.2665e-10;

b1 = c11 + c12.*x + c13.*x.^2 + c14.*x.^3;
b2 = c21 + c22.*x + c23.*x.^2 + c24.*x.^3;
b3 = c31 + c32.*x + c33.*x.^2 + c34.*x.^3;
b4 = c41 + c42.*x + c43.*x.^2 + c44.*x.^3;

h11min = 2.4354e-3;
h11max = 2.5450e-3;
h12min = -6.0339e-5;
h12max = -7.1807e-5;
h13min= 2.1951e-3;
h13max = 1.4580e-3;
h14min = 6.6767;
h14max = 6.9150;
h15min = 9.3228e-1;
h15max = 9.9366e-1;
h21min = 7.7872e-3;
h21max = 7.6828e-3;
h22min = -9.5771e-6;
h22max = -2.4119e-6;
h23min = 6.2229e-4;
h23max = 6.6411e-4;
h24min = 7.7842;
h24max = 7.7461;
h25min = 1.8502;
h25max = 1.9431;
h31min = 9.6340e-1;
h31max = 9.7353e-1;
h32min = 1.5974e-3;
h32max = 1.0577e-3;
h33min = -7.1179e-2;
h33max = -2.1383e-2;
h34min = 2.2320;
h34max = 3.0058;
h35min = 7.8800e-1;
h35max = 9.1845e-1;
h41min = 7.8132e-3;
h41max = 7.3482e-3;
h42min = 9.7085e-11;
h42max = 2.5598e-5;
h43min = 8.2392e-4;
h43max = 1.2457e-3;
h44min = 8.5138;
h44max = 8.1896;
h45min = 2.3125;
h45max = 2.9368;

h51 = 1.9100e-1;
h52 = 7.0300e-2;
h53 = -6.4500e-1;
h54 = 2.0300;
h55 = 1.3000;
h61 = 5.7100e-4;
h62 = 6.1300e-6;
h63 = 5.4700e-4;
h64 = 1.1100;
h65 = 8.3700e-1;

% Combine primary and secondary spectra

for a = 1:length(Rc)
    Elis = Etoa + s(a).*Z./A;
    Beta = sqrt(1-(Ep./(Ep + Elis.*A)).^2); % Particle speed relative to light
    Rlis = 0.001.*sqrt((A.*Elis).^2 + 2.*A.*Ep.*Elis)./Z;
    C = a7 + a8./(1 + exp((Elis - a9)./a10));

    phiTOA = (C.*(Beta.^a5)./(Rlis.^a6)).*(Rtoa./Rlis).^2;
    phiPri = (U./Beta).*phiTOA.*(a2.*exp(-a3.*x) + (1 - a2).*exp(-a4.*x));

    g1min = h11min + h12min.*Rc(a) + h13min./(1 + exp((Rc(a) - h14min)./h15min));
    g1max = h11max + h12max.*Rc(a) + h13max./(1 + exp((Rc(a) - h14max)./h15max));
    g2min = h21min + h22min.*Rc(a) + h23min./(1 + exp((Rc(a) - h24min)./h25min));
    g2max = h21max + h22max.*Rc(a) + h23max./(1 + exp((Rc(a) - h24max)./h25max));
    g3min = h31min + h32min.*Rc(a) + h33min./(1 + exp((Rc(a) - h34min)./h35min));
    g3max = h31max + h32max.*Rc(a) + h33max./(1 + exp((Rc(a) - h34max)./h35max));
    g4min = h41min + h42min.*Rc(a) + h43min./(1 + exp((Rc(a) - h44min)./h45min));
    g4max = h41max + h42max.*Rc(a) + h43max./(1 + exp((Rc(a) - h44max)./h45max));

    phiPmin = g1min.*(exp(-g2min.*x) - g3min.*exp(-g4min.*x)); %length of Rc
    phiPmax = g1max.*(exp(-g2max.*x) - g3max.*exp(-g4max.*x)); %length of Rc

    g5 = h51 + h52.*Rc(a) + h53./(1 + exp((Rc(a) - h54)./h55));
    g6 = h61 + h62.*Rc(a) + h63./(1 + exp((Rc(a) - h64)./h65));

    f3 = g5 + g6.*x;
    f2 = (phiPmin - phiPmax)./(smin.^f3 - smax.^f3);
    f1 = phiPmin - f2.*smin.^f3;

    phiP = f1 + f2.*s(a).^f3;

    phiSec = (phiP.*b1.*E.^b2)./(1 + b3.*E.^b4);

    Ec = (sqrt((1000.*Rc(a).*Z).^2 + Ep.^2) - Ep)./A;
    Es = a13.*(Ec - a14.*x);
    Es1 = max(a15,Es);
    Es2 = max(a16,Es);

    phiPtot = phiPri.*(tanh(a11.*(E./Es1 - 1)) + 1)./2 + ...
        phiSec.*(tanh(a12.*(1 - E./Es2)) + 1)./2;

    clipindex = find(E <= 1, 1, 'last' ); %Make sure the clip index is consistent with the definition of E above

end

P.E = E;
P.phiPtot(i,:) = phiPtot./phiPtot(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot protons spectra and excitation functions
hold on
box on

  [AX1,H1,H2] = plotyy(E,phiPtot./phiPtot(1,1)./P.phiPtot(1,:),E,consts.FepxCl36,@semilogx,@semilogy);
 set(H1,'linewidth',2,'color',[.5,.5,.5])
 set(H2,'linestyle',':','linewidth',2,'color','b','linestyle','--')

  [AX2,H1,H5] = plotyy(E,phiPtot./phiPtot(1,1)./P.phiPtot(1,:),E,consts.O16pxBe10,@semilogx,@semilogx);
 set(H1,'linestyle','-','linewidth',2,'color','k')
 set(H5,'linewidth',2,'color','k','linestyle','-.')

  [AX2,H3,H4] = plotyy(E,phiPtot./phiPtot(1,1)./P.phiPtot(1,:),E,consts.KpxCl36,@semilogx,@semilogx);
 set(AX2(1),'ycolor',[0,0,0],'fontsize',18,'XMinorTick','on','linewidth',3)
 set(AX2(2),'ycolor','k','fontsize',18,'XMinorTick','on','linewidth',3)
 set(H4,'linestyle',':','linewidth',2,'color',[.9,.5,.2])
 set(H3,'linewidth',2,'color',i*[.2,.2,.2])
 set(AX2(1),'YLim',[1 5.5])
 set(AX2(2),'YLim',[.1 1000])
 set(AX2,'XLim',[0 100000])
 xlabel('Energy (MeV)','fontsize',20)
 ylabel(AX1(1),'Normalized Proton Flux')
 ylabel(AX1(2),'millibarns')
 legend([H2,H5,H4,H1],'Fe(p,x)^{36}Cl','^{16}O(p,x)^{10}Be','K(p,x)^{36}Cl','normalized flux','Location','southwest')
 legend boxoff
 text(20000,4.12,'4300 m','fontsize',18)
 text(20000,3.3,'3300 m','fontsize',18)
 text(20000,2.05,'1700 m','fontsize',18)

end

%Print results

set(gcf, "paperunits", "points", "papersize", [900, 1200], 'PaperPosition', [0 0 900 1200]);
print -dpdf -color Normalized_Fluxes.pdf
print -dpng -color Plot_Normalized_Fluxes.png


