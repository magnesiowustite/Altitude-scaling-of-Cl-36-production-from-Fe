% Code for extracting sample parameters into sp
% Modified from Marrero et al. (2016) to include both analytical water
% and pore water

function sp=samppars36_NC(sampledata)
%
% Make sure that sampledata is a vector rather than an array!
%
if (min([size(sampledata,1); size(sampledata,2)]) > 1)
  error('sampledata must be a vector!');
end
%
% Make sampledata a column vectors if it isn't already.
%
if (size(sampledata,1)==1)
  sampledata=sampledata';
end
%
% First, check that the input data is reasonable.
%
if (length(sampledata) ~= 38)
  error('sampledata has wrong size!');
end
%
% Setup the values of sample parameters.
%
sp.concentration36=sampledata(1);
sp.inheritance36=sampledata(2);
sp.epsilon=sampledata(3);
sp.qavg=sampledata(4); %volumetric water/water density
sp.concpore = sp.qavg/(1+sp.qavg); %assume a soil density of 1 for volcaniclastics
% the water content is assumed to be highly uncertain, so knowing the exact
% soil density is not important
sp.LOI=sampledata(5); %water in LOI

%calculate the normalization factor for the oxides after converting the
%fractional water content to weight % water
sp.oxidenormfactor=1-(sp.concpore);

%Corrected water content
sp.qtotal = sp.oxidenormfactor*(sp.LOI/100)+sp.concpore;

sp.ls=sampledata(6);
sp.latitude=sampledata(7);
sp.longitude=sampledata(8);
sp.ST=sampledata(11);
sp.Lambdafe=sampledata(12);
sp.originaloxideinput=sampledata(13:23);

%renormalize the original oxide info to account for water content of rock
newoxidepercent=sp.originaloxideinput*sp.oxidenormfactor;
sp.ci=[NaN; NaN; newoxidepercent;  sampledata(24:31)];

%no renormalization for the target element composition
sp.tci=sampledata(32:36);

% Extract the pressure and elevation from the sampledata vector.
sp.elevation=sampledata(9);
sp.P=sampledata(10);

% Extra: the depth to top of sample - not relevant here
sp.depthtotop=sampledata(37);

sp.tfinal=(sampledata(38)-2010)/1000;
