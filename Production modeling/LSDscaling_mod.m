function [out] = LSDscaling_mod(h,Rc,SPhi,w,consts,nuclide)

% This function is modified from Lifton et al. (2014) to
% calculate reaction-specific scaling factors for all major Cl-36 production
% reactions and Be-10 in quartz.

%-------------------------------------------------------------------------------
% Implements the Lifton Sato Dunai scaling scheme for spallation.
%
% Syntax: scalingfactor = LSDscaling_mod(h,Rc,SPhi,w,consts);
%
% Where:
%   h = atmospheric pressure (hPa)
%   Rc = cutoff rigidity (GV)
%   SPhi = solar modulation potntial (Phi, see source paper)
%   w = fractional water content of ground (nondimensional)

mfluxRef = consts.mfluxRef;
muRef = (mfluxRef.neg + mfluxRef.pos);

% Select reference values for nuclide of interest or flux

if nuclide == 3
    HeRef = consts.P3nRef + consts.P3pRef;
elseif nuclide == 10
    BeRef = consts.P10nRef + consts.P10pRef;
elseif nuclide == 14
    CRef = consts.P14nRef + consts.P14pRef;
elseif nuclide == 26
    AlRef = consts.P26nRef + consts.P26pRef;
else
    SpRef = consts.nfluxRef + consts.pfluxRef;
    % Sato et al. (2008) Reference hadron flux integral >1 MeV
end

EthRef = consts.ethfluxRef;
ThRef = consts.thfluxRef;

% Site nucleon fluxes

NSite = Neutrons_mod(h,Rc,SPhi,w,consts,nuclide);
%[ethflux,thflux] = NeutronsLowE(h,Rc,SPhi,w);
PSite = Protons_mod(h,Rc,SPhi,consts,nuclide);

load ref.mat;
load Reference.mat;

%Nuclide-specific scaling factors as f(Rc)
if nuclide == 10
    out = (NSite.P10n + PSite.P10p)./BeRef;
elseif nuclide == 56
    out = (NSite.PFe36n + PSite.PFe36p)./ref.mat.Fe;
elseif nuclide == 48
    out = (NSite.PTi36n + PSite.PTi36p)./ref.mat.Ti;
elseif nuclide == 39
    out = (NSite.PK36n + PSite.PK36p)./ref.mat.K;
elseif nuclide == 40
    out = (NSite.PCa36n + PSite.PCa36p)./ref.mat.Ca;
elseif nuclide = ('eth')
    %out = ethflux./EthRef; %Epithermal neutron flux scaling factor as f(Rc)
    out = ((NSite.nflux + PSite.pflux))./SpRef;
elseif nuclide = ('th')
    %out = thflux./ThRef; %Thermal neutron flux scaling factor as f(Rc)
    %else    %Total nucleon flux scaling factors as f(Rc)
    out = ((NSite.nflux + PSite.pflux))./SpRef; % Sato et al. (2008) Reference hadron flux integral >1 MeV
end


