function muon = Muon_model(depths,h,Rc,s)

% Model for calculating site-specific single-exponential-term approximation for muons, fit to
% Heisinger model to ~5000 g/cm^2 in the subsurface.

% takes following inputs:
% h = pressure
% Rc = cutoff rigidity
% s = solar modulation constant

pp = physpars_NC; %load physical parameters
z = depths; %setup linespace to 10,000 g/cm^2

muon = muonfluxsato(z,h,Rc,s,pp,'yes'); %calculate muons fluxes
Pmun=pp.Ys*(muon.R)+0.0000058*(muon.phi); %calculate muon-produced neutrons

%Data and model for parameterization
X = z'; %setup fitting depths
P_mu = Pmun';

leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2)); %function to fit
F = leasqrfunc;

%Generate Fit
pin = [10,1500]; % initial value estimates
[f,p,r21]  = leasqr(X,P_mu,pin,F); %least-squares fit

%Outputs
muon.a = p(1,1); %pre-factor
muon.b = p(2,1); %attenuation length
r21;

%Graph model
## close all
## hold all
## hx = loglog(Pmun,X,'linestyle','--','color','b','linewidth',3);
## hy = loglog(f,X,'color','r','linewidth',3,'linestyle','-');
## set (gca (), "ydir", "reverse")
