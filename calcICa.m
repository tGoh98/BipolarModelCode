% Calcium current
function [t, ICa] = calcICa(tspan, tV, voltage, Cas_values)
    AmCa = @(t) 12000*(120-getV(tV,voltage,t))/(exp((120-getV(tV,voltage,t))/25)-1);
    BmCa = @(t) 40000/(exp((getV(tV,voltage,t)+68)/25)+1);
    [t,p] = ode45(@(t,y) AmCa(t)*(1-y)-BmCa(t)*y,tspan,0.13);
    Cao = 2500;
    gca_ = 1.1;
    hCa = @(t) exp(-(getV(tV,voltage,t)-50)/11)/(exp(-(getV(tV,voltage,t)-50)/11)+1);
    hCaVal = zeros(size(t));
    for i = 1:length(t)
         hCaVal(i) = hCa(t(i));
    end
    gCa = gca_*p.^4.*hCaVal;
    ICa = zeros(1, length(t));
    for i = 1:length(t)
        Cas = getCas(tV,Cas_values,i);
%         Cas = 2500;
        ECa = 12.9*log(Cao/Cas);
        ICa(i) = gCa(i)*(getV(tV,voltage,i)-ECa);
    end
end