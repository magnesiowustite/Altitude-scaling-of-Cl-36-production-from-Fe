function out = air_pressure(lat,lon,alt,atm)

% Calculates air pressure at an input location.
% Modified from Lifton et al. (2014).

%-------------------------------------------------------------------------------
% syntax : LSD_mod(lat,lon,alt,atm);

% Where:
% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes,
%     or 0-360 degrees E)
% alt = altitude (m asl)
% atm = 1 for standard atmosphere, atm = 2 for ERA40 model

% Load the input data structure
sample.lat = lat;
sample.lon = lon;
sample.alt = alt;
sample.atm = atm;

%Determine pressure at elevation
if sample.atm == 1
    stdatm = 1;
    gmr = -0.03417; % Assorted constants
    dtdz = 0.0065; % Lapse rate from standard atmosphere
else
    stdatm = 0;
end

% Pressure correction
if stdatm == 1
    % Calculate site pressure using the Standard Atmosphere parameters with the
    % standard atmosphere equation.
    sample.pressure = 1013.25 .* exp( (gmr./dtdz) .* ( log(288.15) - log(288.15 - (alt.*dtdz)) ) );
else
    sample.pressure = ERA40atm(sample.lat,sample.lon,sample.alt);
end

%output pressure
out = sample.pressure;



