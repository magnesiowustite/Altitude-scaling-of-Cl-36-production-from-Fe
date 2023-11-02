%Function to fit an exponential approximation to the polynomials for Cl-36 from
%K and Be-10 in quartz attenuation in granite from Argento (2015).

function Ratio = Argento_attenuation()

%Setup depth profile
z = linspace(0,500);
%Calculate production rates with depth
P_K = (0.000000008834*z.^3-0.000005795*z.^2-0.005548*z+4.954); %Cl-36 from K
P_K = e.^P_K;
P_Be = (0.0000000005952*z.^3-0.0000008528*z.^2-0.006895*z+1.52); %Be-10 in qtz
P_Be = e.^P_Be;

%Plot ratio
##hold on
##plot((P_K./P_Be)./(P_K(1,1)./P_Be(1,1)),z)
##set (gca (), "ydir", "reverse")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data and model for parameterization for K
X = z'; %setup fitting depths
P_K = P_K';

leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2)); %function to fit
F = leasqrfunc;

%Generate Fit
pin = [100,160]; % initial value estimates
[f,p,r21]  = leasqr(X,P_K,pin,F); %least-squares fit

%Outputs
a_K = p(1,1); %pre-factor
b_K = p(2,1); %attention length
r21;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data and model for parameterization for Be
X = z'; %setup fitting depths
P_K = P_Be';

leasqrfunc = @(x,p,r21) p(1)*exp(-X/p(2)); %function to fit
F = leasqrfunc;

%Generate Fit
pin = [100,160]; % initial value estimates
[f,p,r21]  = leasqr(X,P_K,pin,F); %least-squares fit

%Outputs
a_Be = p(1,1); %pre-factor
b_Be = p(2,1); %attention length
r21;

%Calculate attenuation length ratio
Ratio = b_K./b_Be;

%Plot best-fit model
##P_K = a_K*e.^(-z./b_K);
##P_Be = a_Be*e.^(-z./b_Be);
##plot((P_K./P_Be)./(P_K(1,1)./P_Be(1,1)),z);



