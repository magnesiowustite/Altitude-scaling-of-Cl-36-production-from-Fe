function out = avg_RC(lat,lon,age)

% This function makes linear RC and SPhi vectors extending
% from the present to the input age. Modified from Lifton et al. (2014).

%------------------------------------------------------------------------------
% syntax : avg_RC(lat,lon,age);

% Where:
% lat = sample latitude in deg N (negative values for S hemisphere)
% lon = sample longitude in deg E (negative values for W longitudes,
%     or 0-360 degrees E)
% age = endpoint (years BP)

load consts_LSD.mat;

% Load the input data structure
sample.lat = lat;
sample.lon = lon;
sample.age = age;

% Make the time vector

% Age Relative to t0=2010
tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];
LSDRc = zeros(1,length(tv));

% Need solar modulation parameter
this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

% interpolate an M for tv > 7000...
temp_M = interp1(consts.t_M,consts.M,tv(77:end));

% catch for negative longitudes before Rc interpolation
if sample.lon < 0; sample.lon = sample.lon + 360;
  end

% Make up the Rc vectors.
[loni,lati,tvi] = meshgrid(sample.lon,sample.lat,tv(1:76));
LSDRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,loni,lati,tvi);

% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
dd = [6.89901,-103.241,522.061,-1152.15,1189.18,-448.004;];

LSDRc(77:end) = temp_M.*(dd(1)*cosd(sample.lat) + ...
   dd(2)*(cosd(sample.lat)).^2 + ...
   dd(3)*(cosd(sample.lat)).^3 + ...
   dd(4)*(cosd(sample.lat)).^4 + ...
   dd(5)*(cosd(sample.lat)).^5 + ...
   dd(6)*(cosd(sample.lat)).^6);

% Next, chop off tv
clipindex = find(tv <= sample.age, 1, 'last' );
tv2 = tv(1:clipindex);
if tv2(end) < sample.age;
    tv2 = [tv2 sample.age];
end;

% Now shorten the Rc's commensurately
LSDRc = interp1(tv,LSDRc,tv2);
LSDSPhi = interp1(tv,this_SPhi,tv2);

%interpolate Rc and SPhi onto a uniform linear scale
out.LSDRc = interp1(tv2,LSDRc,linspace(0,sample.age,sample.age/100));
out.LSDSPhi = interp1(tv2,LSDSPhi,linspace(0,sample.age,sample.age/100));



