function []=Nuclide_Altitude_Plot()

%Makes figure 2, a plot of production ratios with altitude at high latitude
clf

load consts_LSD.mat;

% Load the input data structure

for i = 490:10:1013.25

sample.lat = 90;
sample.lon = 0;
sample.pressure = i;
sample.age = 0;
w = 0.066; %average gravimetric water content of the ground

% Make the time vector
calFlag = 0;

% Age Relative to t0=2010
tv = [0:10:50 60:100:50060 51060:1000:2000060 logspace(log10(2001060),7,200)];

LSDRc = zeros(1,length(tv));

% Need solar modulation parameter
this_SPhi = zeros(size(tv)) + consts.SPhiInf; % Solar modulation potential for Sato et al. (2008)
this_SPhi(1:120) = consts.SPhi; % Solar modulation potential for Sato et al. (2008)

if w < 0
    w = 0.066; % default gravimetric water content for Sato et al. (2008)
end

% interpolate an M for tv > 7000...
temp_M = interp1(consts.t_M,consts.M,tv(77:end));

% catch for negative longitudes before Rc interpolation
if sample.lon < 0; sample.lon = sample.lon + 360; end;

% Make up the Rc vectors.

% Modified to work with new interpolation routines in MATLAB 2012a and later. 09/12
[loni,lati,tvi] = meshgrid(sample.lon,sample.lat,tv(1:76));
LSDRc(1:76) = interp3(consts.lon_Rc,consts.lat_Rc,consts.t_Rc,consts.TTRc,loni,lati,tvi);

% Fit to Trajectory-traced GAD dipole field as f(M/M0), as long-term average.
%
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

LSDout_56 = LSDscaling_mod(sample.pressure,LSDRc(:),LSDSPhi,w,consts,56);
LSDout_39 = LSDscaling_mod(sample.pressure,LSDRc(:),LSDSPhi,w,consts,39);
LSDout_10 = LSDscaling_mod(sample.pressure,LSDRc(:),LSDSPhi,w,consts,10);

PRR_56 = mean(LSDout_56)/mean(LSDout_10);
PRR_39 = mean(LSDout_39)/mean(LSDout_10);
PRR_10 = mean(LSDout_10)/mean(LSDout_10);

out_56(1,54-(i-480)/10) = PRR_56;
out_39(1,54-(i-480)/10) = PRR_39;
out_10(1,54-(i-480)/10) = PRR_10;

end

out_56 = out_56./out_56(1,1);
out_39 = out_39./out_39(1,1);
out_10 = out_10./out_10(1,1);

%Normalize Output to SLHL

pressure = [1013.25:-10:490];

%Compute Elevations

for j = 1:53

elevation(1,j) = (288.15-288.15*(pressure(1,j)/1013.25)^(1/(.03417/.0065)))/.0065/1000;

end

% Make plot

hold all
H1 = plot(out_56,pressure,'color',[0,0,.75],'linewidth',2);
set(gca,'Ydir','reverse',"fontsize", 18,'fontweight','normal', "linewidth", 3, 'XMinorTick', 'on','YMinorTick','on','box','on');
H2 = plot(out_39,pressure,'color',[.8,.2,.8],'linewidth',2);
H3= plot(out_10,pressure,'color','k','linewidth',2);
H4 = plot(out_56./out_39,pressure,'color',[.9,.5,.2],'linewidth',2);
plot([0,2],[602.97,602.97],'linewidth',2,'linestyle','--','color','k')
plot([0,2],[683.19,683.19],'linewidth',2,'linestyle','--','color','k')
plot([0,2],[830.94,830.94],'linewidth',2,'linestyle','--','color','k')

text(.6,590,'Mt. Evans - 4300 m','fontsize',18,'fontweight','normal')
text(.6,670,'Sierra Nevada - 3300 m','fontsize',18,'fontweight','normal')
text(.6,817,'Owens Valley - 1700 m','fontsize',18,'fontweight','normal')

legend([H2,H3,H1,H4],'^{36}Cl_{K}/^{10}Be_{qtz}','^{10}Be_{qtz}/^{10}Be_{qtz}','^{36}Cl_{Fe}/^{10}Be_{qtz}','^{36}Cl_{Fe}/^{36}Cl_{K}','^{36}Cl/^{10}Be','Location','southwest')
legend boxoff

hold off

xlim([.5,1.5]);
ylim([500,1013.25]);

xlabel('Production ratio deviation from sea level',"fontweight","normal","fontsize",20);
ylabel ('Pressure (hPa)',"fontweight","normal","fontsize",20);

% Write results to file
set(gcf, "paperunits", "points", "papersize", [900, 600], 'PaperPosition', [0 0 900 600]);
print -dpdf -color altitude_model.pdf
print -dpng -color altitude_model.png
