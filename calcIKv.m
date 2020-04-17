% Delayed rectifying potassium current
function [t, IKv] = calcIKv(tspan, tV, voltage)
    [t,p] = ode45(@(t,y) calcProb(t,y,tV,voltage),tspan,[0.5,0]);
    gkv_ = 2.0;
    Ek = -58;
    gKv = gkv_*p(:,1).^3.*p(:,2);
    IKv = zeros(1, length(t));
    for i = 1:length(t)
        IKv(i) = gKv(i)*(getV(tV,voltage,t(i))-Ek);
    end
end

function probabilities = calcProb(t,y,tV,voltage)
    AmKv = @(t) 400/(exp(-(getV(tV,voltage,t)-15)/36)+1);
    BmKv = @(t) exp(-getV(tV,voltage,t)/13);
    AhKv = @(t) 0.0003*exp(-getV(tV,voltage,t)/7);
    BhKv = @(t) 80/(exp(-(getV(tV,voltage,t)+115)/15)+1)+0.02;
    probabilities = zeros(2,1);
    probabilities(1) = AmKv(t)*(1-y(1))-BmKv(t)*y(1);
    probabilities(2) = AhKv(t)*(1-y(2))-BhKv(t)*y(2);
end