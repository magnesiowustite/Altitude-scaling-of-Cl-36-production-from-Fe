%% Function for setting-up muon parameters

function [out] = muon_parameters()
  %load input data
  sampledata = xlsread('Mt Evans.xlsx','Inputs','B2:AM49');

  for i = 1:rows(sampledata)

  %setup data
  data    = sampledata(i,:); %select data for sample

  %calculate denudation rate
  pp=physpars_NC();
  sp36 = samppars36_NC(data);
  cp36 = comppars36_muons(pp,sp36);
  cp10 = comppars1026_muons(pp,sp36);

  %put production fractions in output
  out.P_neg36_a = cp36.P_negtotal36_a; %Cl-36 slow production rate a
  out.P_neg36_c = cp36.P_negtotal36_c; %Cl-36 slow production rate c
  out.P_fast36 = cp36.P_fasttotal36; %Cl-36 fast production rate
  out.P_neg10_a = cp10.P_neg10_a; %Be-10 slow production rate a
  out.P_neg10_c = cp10.P_neg10_c; %Be-10 slow production rate c
  out.P_fast10 = cp10.P_fast10; %Be-10 fast production rate
  out.b = cp36.b; %slow attenuation 1
  out.d = cp36.d; %slow attenuation 2
  out.f = cp36.f; %fast attenuation

  end
