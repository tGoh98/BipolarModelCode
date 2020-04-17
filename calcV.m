function [t, res] = calcV(tspan)
[t, res] = ode45(@(t,y) solveEqs(t,y),tspan,[0.5,0,0.13,0.37,-30,0,0,0,0,0,2500,2500]);
end

function results = solveEqs(t,y)
% 12 total differential equations
results = zeros(12,1);
C = 0.01;

% Intermediate values
mKv = y(1);
hKv = y(2);
mCa = y(3);
mKc = y(4);
V = y(5);
C1 = y(6);
C2 = y(7);
O1 = y(8);
O2 = y(9);
O3 = y(10);
Cas = y(11);
Cad = y(12);

% Delayed rectifying potassium current
gkv_ = 2.0;
Ek = -58;
AmKv = 400/(exp(-(V-15)/36)+1);
BmKv = exp(-V/13);
AhKv = 0.0003*exp(-V/7);
BhKv = 80/(exp(-(V+115)/15)+1)+0.02;
gKv = gkv_*mKv^3*hKv;
IKv = gKv*(V-Ek);

% Hyperpolarization activated current
gh_ = 0.975;
Eh = -17.7;
Ah = 3/(exp((V+110)/15)+1);
Bh = 1.5/(exp(-(V+115)/15)+1);
mh = O1+O2+O3;
gh = gh_*mh;
Ih = gh*(V-Eh);

% Calcium current
Cao = 2500;
gca_ = 1.1;
AmCa = 12000*(120-V)/(exp((120-V)/25)-1);
BmCa = 40000/(exp((V+68)/25)+1);
ECa = 12.9*log(Cao/Cas);
hCa = exp(-(V-50)/11)/(exp(-(V-50)/11)+1);
gCa = gca_*mCa^4*hCa;
ICa = gCa*(V-ECa);

% Ca-dependent K Current
gkc_ = 8.5;
Ek = -58;
AmKc = 100*(230-V)/(exp((230-V)/52)-1);
BmKc = 120*exp(-V/95);
mKc1 = Cas/(Cas+0.2);
gKc = gkc_*mKc^2*mKc1;
IKCa = gKc*(V-Ek);

% Leakage current
gl = 0.23;
El = -21;
Il = gl*(V-El);

% Constants/parameters for Cas and Cad
F = 9.649e5;
Vs = 1.692e-13;
Vd = 7.356e-13;
Dca = 6e-8;
Ssd = 4e-8;
dsd = 5.9e-5;

% Results
results(1) = AmKv*(1-mKv)-BmKv*mKv; %mKv
results(2) = AhKv*(1-hKv)-BhKv*hKv; %hKv
results(3) = AmCa*(1-mCa)-BmCa*mCa; %mCa
results(4) = AmKc*(1-mKc)-BmKc*mKc; %mKc
results(5) = (-(IKv+Ih+ICa+IKCa+Il)+getI(t))/C; %V
results(6) = -C1*4*Ah+C2*Bh; %C1
results(7) = -C1*4*Ah-C2*(3*Ah+Bh)+O1*2*Bh; %C2
results(8) = C2*3*Ah-O1*(2*Ah+2*Bh)+O2*3*Bh; %O1
results(9) = O1*2*Ah-O2*(Ah+3*Bh)+O3*4*Bh; %O2
results(10) = O2*Ah-O3*4*Bh; %O3
results(11) = -ICa/(2*F*Vs)-((Dca*Ssd)/(Vs*dsd))*(Cas-Cad); %Cas || ignoring calcium pump and exchanger for now
results(12) = ((Dca*Ssd)/(Vd*dsd))*(Cas-Cad); %Cad
end