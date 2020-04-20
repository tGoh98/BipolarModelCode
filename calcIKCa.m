% Calcium-dependent K current
function [t, IKCa] = calcIKCa(tspan, tV, voltage, Cas_values)
    AmKc = @(t) 100*(230-getV(tV,voltage,t))/(exp((230-getV(tV,voltage,t))/52)-1);
    BmKc = @(t) 120*exp(-getV(tV,voltage,t)/95);
    [t,p] = ode45(@(t,y) AmKc(t)*(1-y)-BmKc(t)*y,tspan,0.37);
    gkc_ = 8.5;
    Ek = -58;
    IKCa = zeros(1, length(t));
    for i = 1:length(t)
%         Cas = 2500; % don't know this value
        Cas = getCas(tV,Cas_values,i);
        mKc1 = Cas/(Cas+0.2);
        gKc = gkc_*p(i)^2*mKc1;
        IKCa(i) = gKc*(getV(tV,voltage,i)-Ek);
    end
end