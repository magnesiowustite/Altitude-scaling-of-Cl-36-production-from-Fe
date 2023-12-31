% Code from Marrero et al. (2016) that sets up parameters for Cl-36 dating.
% Only values for the Sa scaling model are included here. Tables were updated...
% to reverse Cr and Li order and include updated parameter values...
% (after Marrero et al., 2021).

function pp=physpars_NC()

% decay constants
pp.lambda36Cl=2.3028e-6;
pp.lambda10Be=4.998e-7;
pp.lambda14C=1.21e-4;
pp.lambda3He=0.0;
pp.lambda21Ne=0.0;
pp.lambda26Al=9.83e-7;

%parameters for chlorine IDMS
pp.Wn=3.54527e+1;
pp.Wo=3.4963e+1;
pp.Sn=3.12e+0;
pp.So=2.93e+2;
pp.An=pp.Sn/(pp.Sn+1);
pp.Ao=pp.So/(pp.So+1);

%
% The table below corresponds to the table on the TABLE page in CHLOE.  The data is from Fabryka-Martin 1988, and was updated by Vasily Alfimov using the
% "Atlas of Neutron Resonances" by S.F. Mughabghab, Elsevier Science, 2006.The rows correspond to the following elements (top to bottom):
%  O, H, C, Na, Mg, Al, Si, P, K, Ca, Ti, Mn, Fe, Cl, B, Sm, Gd, U, Th, Li,
%  Cr.  Li & Cr are from the original Fabryka-Martin table and were not
%  updated.
% The columns correspond to the following categories: col 1= Ai=Atomic weight of element (g/mol),
%  col2 = average log decrement of energy per neutron collison with element i [-], col 3= neutron scattering cross-section of element i[i],
%  col 4= thermal neutron absorption cross-section of element i [-], col 5= dilute resonance integral for element i [-],
%  col 6= Mass stopping power of element i for alpha particles of a given energy [-],
%  col 7= neutron yield of element i per ppm U in radioequilibrium [-], col 8= Neutron yield of element i per ppm Th in radioequilibrium,
%  col 9 = Km, 602/atomic weight of sample, used to convert to at/g from
%  ppm. col 10= Stoichiometry ratio of oxide (oxygen to element i) [-], col 11=atomic number [-],
% col 12=Average capture probability relative to oxygen P(Z) (ref: von
% Egidy & Hartman 1982, value for C from back calculation from Heisinger
% 2002 compound factor calculations).
%
pp.table=[...
16          0.120004104		3.761		0.00019	0.0002693   539     0.23	0.079	0		0   8   1.00;...
1.01		1               20.49		0.3326	0           1542	0       0		0		0.5 1   0.00;...
12.01		0.157760474		4.74		0.0035	0.0018      573     0.45	0.18	13.691	2   6   0.36;...
22.99		0.084543589		3.038		0.517	0.311		454     14.5	6.8		19.42	0.5 11  1.00;...
24.31		0.08009077		3.414		0.0666	0.038		463     5.8     2.6		14.94	1   12  0.93;...
26.98		0.072337427		1.413		0.231	0.17		449     5.1     2.6		11.812	1.5 13  0.76;...
28.09		0.069559975		2.044		0.171	0.082		455     0.69	0.335	10.01	2   14  0.84;...
30.97		0.063210393		3.134		0.165	0.079		444     0       0		8.48	2.5 15  1.04;...
39.1		0.050295528		2.04		2.1		1           432     0.45	0.305	12.78	0.5 19  1.54;...
40.08		0.049086179		2.93		0.43	0.233		436     0       0		10.73	1   20  1.90;...
47.87		0.041208508		4.09		6.41	3.1         367     0       0		7.53	2   22  2.66;...
54.94		0.035968204		2.06		13.36	13.4		351     0       0		8.486	1   25  2.73;...
55.85		0.035390922		11.35		2.56	1.36		353     0.19	0.205	7.54	1.5 26  3.28;...
35.45		0.055371497		15.8		33.14	13.83		420     0       0		16.98	0   17  1.32;...
10.81		0.174236264		4.27		767		343         537     62.3	19.2	55.68	0   5   0.25;...
150.36      0.013242694		38  		9640	1400		0       0       0		4.004	0   62  4.4;...
157.25      0.012664667		172     	41560	390         0       0       0		3.828	0   64  5.8;...
238.03      0.008378873		9.08		2.68	277         0       0       0		2.529	0   92  4.7;...
232.04      0.008594581		13.55		7.34	83.3		0       0       0		2.594	0   90  3.0;...
52.0        0.038           3.38        3.1     1.6         0.0     0.0     0.0     11.578  0.0 24  2.98;...
6.9         0.264           0.95        70.5    0.0         548     21.1    9.6     86.731  0.0 3   0.18];
%reversed Cr and Li order as in the corrigendum to Marrero et al. 2016. 8.18.22

pp.ag=6.02e23;             % Avagadro's number.
pp.Sigmaetha=0.0548;
pp.Sigmatha=0.060241;
pp.Sigmasca=0.3773;
pp.Aa=14.5;
%
% Published values were used for these parameters.
%
pp.PsTi0=6.5002E-22;
%pp.PsTi0=3.98002E-22*3.5; Value from Fink (2000)
pp.PsFe0=1.28e-22; %value from Stone (2005) abstract
pp.YstarCa=2.0e-24;
pp.YstarK=3.85e-24;
pp.Psimu0=175.0;

pp.Natoms3 = 2.006e22;
pp.Natoms10 = 2.006e22;
pp.Natoms14 = 2.006e22;
pp.Natoms26 = 1.003e22;

%k_negpartial must be multiplied
%by fstar(below)
pp.k_negpartial14 = 0.704 * 0.1828;
pp.delk_negpartial14 = 0.704 * 0.1828 * 0.0011;
pp.k_negpartial36K=0.802;
pp.k_negpartial36Ca=0.8486;
pp.k_negpartial10=(0.704 * 0.1828)./1.106;
pp.sigmak_neg10 = (0.704 * 0.1828 * 0.0003)./1.106;
pp.k_negpartial26=0.296 * 0.6559;
pp.sigmak_neg26 = 0.296 * 0.6559 * 0.002;

pp.sigma190_14 = 0.45e-27;
pp.delsigma190_14 = 0.25e-27;
pp.sigma190_10 = (0.094e-27)./1.106;
pp.sigmasigma190_10 = (0.013e-27)./1.106;
pp.sigma190_26 = 1.41e-27;
pp.sigmasigma190_26 = 0.17e-27;
pp.sigma190_36Ca=1.40e-27;
pp.sigmasigma190_36Ca=0.30e-27;
%
% More parameters.
%
pp.pEtha=0.56;
pp.La=3.9;
pp.Xa=0.1346;
pp.Phimuf0=700000;
pp.Ys=0.44;
pp.Yf=1.0;

%This value needed for muons and is from Nat's geomag scaling model values
pp.SPhiInf = 416.492; % Changed 12/13/11 to reflect updated SPhi values from Usoskin et al 2011

%
% Calibrated production rates.  See calibratexx.out.
%
%Production rates for all nuclides/all scaling schemes. Uncomment the one
%for the appropriate scaling. Muon parameters calibrated for Al, Be, and
%Cl.

% SF scaling
	pp.PsAl=28.6;%CRONUScalc value - Borchers et al., 2015
	pp.sigmaPsAl=3.3;%CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.PsBe=4.09;%CRONUScalc value - Borchers et al., 2015
	pp.sigmaPsBe=0.35;%CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.PsCa0=56.24897;  %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigmaPsCa0=4.6; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.PsK0=153.296853; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigmaPsK0=12; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.Pf0=601; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigmapf0=179; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.PsC=12.7;%CRONUScalc value - Borchers et al., 2015
	pp.sigmaPsC=0.0;%CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.PsHe=119;%CRONUScalc value - Borchers et al., 2015
	pp.sigmaPsHe = 19; %CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
    pp.PsNe=16.96;% From Balco & Shuster 2009 (ratio of Be/Ne)
	pp.sigmaPsNe =2.12; %Value from Balco's paper, error propagated.
	pp.fstar14=0.137; %Not calibrated as part of CRONUS
	pp.fstar10=1.89e-3; %CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.fstar26=11.7e-3; %CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.fstar36K=0.05736245; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.fstar36Ca=0.01358239; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigma014=8.79321E-30; %Not calibrated as part of CRONUS
	pp.sigma036K=9.398043e-30; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigma036Ca=8.261573e-30; %CRONUScalc value - Marrero et al., 2015 (36Cl calibration)
	pp.sigma010=0.252e-30; %CRONUScalc value - Phillips et al., 2015 (Synthesis paper)
	pp.sigma026=4.10e-30; %CRONUScalc value - Phillips et al., 2015 (Synthesis paper)

