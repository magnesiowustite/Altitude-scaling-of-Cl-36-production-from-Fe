function []=Final_Plot()

%% Makes Figure 4, a plot of calibrated scaling factor ratios vs. model predictions
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First calculate the model
load consts_LSD.mat;
w = 0.066; %gravimetric water content of the ground
Rc = 6.3; %GV

% Load the input data structure

for i = 490:10:1013.25

sample.pressure = i;

LSDout_56 = LSDscaling_mod(sample.pressure,Rc,462,w,consts,56);
LSDout_39 = LSDscaling_mod(sample.pressure,Rc,462,w,consts,39);
LSDout_10 = LSDscaling_mod(sample.pressure,Rc,462,w,consts,10);

PRR_56Be = mean(LSDout_56)/mean(LSDout_10);
PRR_56Cl = mean(LSDout_56)/mean(LSDout_39);

out_10(1,54-(i-480)/10) = PRR_56Be;
out_39(1,54-(i-480)/10) = PRR_56Cl;

end
out_10(1,1);
out_39(1,1);

out_10 = out_10./out_10(1,1);
out_39 = out_39./out_39(1,1);

%Normalize Output to SLHL

pressure = [1013.25:-10:490];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now plot the model
hold all
H0 = plot(out_10,pressure./100,'color',[0,0,.75],'linewidth',2);
H1 = plot(out_39,pressure./100,'color',[.9,.5,.2],'linewidth',2);
set(gca,'Ydir','reverse',"fontsize", 18,'fontweight','normal', "linewidth", 3,...
 'XMinorTick', 'on','YMinorTick','on','box','on');

xlim([1,1.5]);
ylim([5,10.1325]);

xlabel('Scaling factor ratio deviation from sea level',"fontweight","normal","fontsize",20);
ylabel ('Pressure (hPa)',"fontweight","normal","fontsize",20);

%Now plot measurement data

%load data
sampledata = xlsread('Input data.xlsx','Results','G5:H25');

%Owens Valley
data_1700_Cl = [sampledata(1,1),sampledata(1,2);
sampledata(2,1),sampledata(2,2)];
data_1700_Be = [sampledata(4,1),sampledata(4,2);
sampledata(5,1),sampledata(5,2)];
%Mt. Evans
data_4300_Cl = [sampledata(8,1),sampledata(8,2)];
data_4300_Be = [sampledata(9,1),sampledata(9,2)];
%Sierra Nevada
data_3300 = [sampledata(12,1),sampledata(12,2);
sampledata(13,1),sampledata(13,2);
sampledata(14,1),sampledata(14,2);
sampledata(15,1),sampledata(15,2);
sampledata(16,1),sampledata(16,2)];
%Averages
Avg_1700_Cl = [sampledata(19,1),sampledata(19,2)];
Avg_1700_Be = [sampledata(20,1),sampledata(20,2)];
Avg_3300 = [sampledata(21,1),sampledata(21,2)];

%pressures
altitudes_1700 = [830.844 830.8438]./100;
altitudes_3300 = [677.5011 677.671 677.9248 680.6452 681.5138]./100;
altitudes_4300 = [602.97]./100;

H4 = plot(data_1700_Cl(:,1),altitudes_1700);
H5 = plot(data_1700_Be(:,1),altitudes_1700);
H6 = plot(data_3300(:,1),altitudes_3300);

set(H4,'color','k','linestyle','none','marker','d','markersize',13,'markerfacecolor',[.9,.9,.9])
set(H5,'color','k','linestyle','none','marker','o','markersize',13,'markerfacecolor',[.9,.9,.9])
set(H6,'color','k','linestyle','none','marker','o','markersize',13,'markerfacecolor',[.9,.9,.9])

%Plot Averages
H9 = errorbar(Avg_1700_Cl(:,1),mean(altitudes_1700),Avg_1700_Cl(:,2),'>');
H10 = errorbar(Avg_1700_Be(:,1),mean(altitudes_1700),Avg_1700_Be(:,2),'>');
H11 = errorbar(Avg_3300(:,1),mean(altitudes_3300),Avg_3300(:,2),'>');
H12 = errorbar(data_4300_Cl(:,1),[6.0297],data_4300_Cl(:,2),'>');
H13 = errorbar(data_4300_Be(:,1),[6.0297],data_4300_Be(:,2),'>');

set(H9,'color','k','linestyle','none','marker','d','markersize',15,'markerfacecolor',[.9,.5,.2],'linewidth',2);
set(H10,'color','k','linestyle','none','marker','o','markersize',12,'markerfacecolor',[0,0,.75],'linewidth',2);
set(H11,'color','k','linestyle','none','marker','o','markersize',12,'markerfacecolor',[0,0,.75],'linewidth',2);
set(H12,'color','k','linestyle','none','marker','d','markersize',15,'markerfacecolor',[.9,.5,.2],'linewidth',2);
set(H13,'color','k','linestyle','none','marker','o','markersize',12,'markerfacecolor',[0,0,.75],'linewidth',2);

set(gca, 'ytick', [10,9,8,7,6,5]);
set(gca, 'yticklabel', { "1000", "900", "800", "700", "600", "500" });

%%%%%%%%%%%%%%%%%%%%%% Fit polynomial coefficients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data and model for parameterization for Be
X = 1013.25-pressure(1:51)'; %setup fitting depths
Ratio = out_10(1:51)';

coeff = polyfit(X,Ratio,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

expo = coeff(1,3)+coeff(1,2)*(1013.25-pressure)+coeff(1,1)*(1013.25-pressure).^2;

H11 = plot(expo,pressure./100,'linewidth',4,'linestyle',':','color','k','marker','none');

legend([H11,H0,H1,H10,H9,H5,H4],'polynomial fit','^{36}Cl_{K}/^{10}Be_{qtz}','^{36}Cl_{Fe}/^{36}Cl_{K}',...
'Site mean, Fe/Be','Site mean, Fe/K','Qtz. Sample','Fs. Sample','Location','southeast','fontsize',12);
legend boxoff

%Write results to file
set(gcf, "paperunits", "points", "papersize", [900, 600], 'PaperPosition', [0 0 900 600]);
%print -dpdf -color Plot_Data_Final.pdf
print -dpng -color Plot_Data_Final.png
