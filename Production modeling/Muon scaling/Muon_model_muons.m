function muon = Muon_model_muons(depths,h,Rc,s)

% Model for calculating site-specific 3-term approximation for muons, fit to
% Heisinger model to ~5000 g/cm^2 in the subsurface.

% takes following inputs:
% h = pressure
% Rc = cutoff rigidity
% s = solar modulation constant

depths = [linspace(0,5000,30)]; %we will fit to 5000 g/cm^2

pp = physpars_NC; %load physical parameters
z = depths;

muon = muonfluxsato(z,h,Rc,s,pp,'yes'); %calculate muons fluxes

%Data and model for parameterization
X = z'; %setup fitting depths
P_mu_R = muon.R'; %negative muons
P_mu_phi = muon.phi'; %fast muons

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fist, fit the negative muons - 2 terms
leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2))+p(3)*exp(-X/p(4)); %function to fit
F = leasqrfunc;

%Generate Fit
 stol=.0001; niter=500;
 weights = ones(size(P_mu_phi));
pin = [10,10,10,10]; % initial value estimates
options.bounds = [0,10000; 0,15000; 0,15000; 0,15000];
[f,p,r21]  = leasqr(X,P_mu_R,pin,F,stol,niter,weights,.001*ones(size(pin)),...
"dfdp",options); %least-squares fit

%Output coefficients and attenuation lengths for negative muons
muon.a = p(1,1);
muon.b = p(2,1);
muon.c = p(3,1);
muon.d = p(4,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Next, fit the fast muons
leasqrfunc = @(x,q,r21) q(1)*exp(-X/q(2)); %function to fit
F = leasqrfunc;

%Generate Fit
 stol=.0001; niter=500;
 weights = ones(size(P_mu_phi));
pin = [10,1500]; % initial value estimates
options.bounds = [0,10E6; 0,10000];
[f,q,r21]  = leasqr(X,P_mu_phi,pin,F,stol,niter,weights,.001*ones(size(pin)),...
"dfdp",options); %least-squares fit

%Output coefficients and attenuation lengths for negative muons
muon.e = q(1,1);
muon.f = q(2,1);

##hold on
##%X
##plot(X,P_mu_phi)
##%y = p(1)*exp(-X/p(2))+p(3)*exp(-X/p(4));
##y = q(1)*exp(-X/q(2));
##plot(X,y,'linestyle','--')
