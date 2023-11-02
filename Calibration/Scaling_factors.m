function out = Scaling_factors()

% This code calculates time-averaged, nuclide-specific scaling factors for
% Be-10 and all Cl-36 production pathways for the locations and exposure ages
% specified in "Input data.xlsx".

inputs = xlsread('Input data.xlsx','Inputs','B14:G22');

for i = 1:rows(inputs)

    lat = inputs(i,1);
    lon = inputs(i,2);
    alt = inputs(i,3);
    atm = inputs(i,4);
    w = inputs(i,5);
    age = inputs(i,6);

    RCs = avg_RC(lat,lon,age);
    pressure = air_pressure(lat,lon,alt,atm);
    load consts_LSD.mat;

  out(i,2) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,56));
  out(i,1) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,10));
  out(i,3) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,48));
  out(i,4) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,40));
  out(i,5) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,39));
  out(i,6) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,'eth'));
  out(i,7) = mean(LSDscaling_mod(pressure,RCs.LSDRc,RCs.LSDSPhi,w,consts,'th'));
  %out(i,8) = ERA40atm(lat,lon,alt);

end

