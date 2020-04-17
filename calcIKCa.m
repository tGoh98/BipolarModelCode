% Calcium-dependent K current
function [t, IKCa] = calcIKCa(tspan, tV, voltage)
    AmKc = @(t) 100*(230-getV(tV,voltage,t))/(exp((230-getV(tV,voltage,t))/52)-1);
    BmKc = @(t) 120*exp(-getV(tV,voltage,t)/95);
    [t,p] = ode45(@(t,y) AmKc(t)*(1-y)-BmKc(t)*y,tspan,0.37);
    Cas = 2500; % don't know this value
    mKc1 = Cas/(Cas+0.2);
    gkc_ = 8.5;
    Ek = -58;
    gKc = gkc_*p.^2*mKc1;
    IKCa = zeros(1, length(t));
    for i = 1:length(t)
        IKCa(i) = gKc(i)*(getV(tV,voltage,i)-Ek);
    end
end