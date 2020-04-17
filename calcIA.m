% Transient outward current
function [t, IA] = calcIA(tspan, tV, voltage)
    [t,p] = ode45(@(t,y) calcProb(t,y,tV,voltage),tspan,[0,0]);
    ga_ = 35;
    Ek = -58;
    gA = ga_*p(:,1).^3.*p(:,2);
    IA = zeros(1, length(t));
    for i = 1:length(t)
        IA(i) = gA(i)*(getV(tV,voltage,t(i))-Ek);
    end
end

function probabilities = calcProb(t,y, tV, voltage)
    AmA = @(t) 1200/(exp(-(getV(tV,voltage,t)-50)/28)+1);
    BmA = @(t) 6*exp(-getV(tV,voltage,t)/10);
    AhA = @(t) 0.045*exp(-getV(tV,voltage,t)/13);
    BhA = @(t) 75/(exp(-(getV(tV,voltage,t)+50)/15)+1);
    probabilities = zeros(2,1);
    probabilities(1) = AmA(t)*(1-y(1))-BmA(t)*y(1);
    probabilities(2) = AhA(t)*(1-y(2))-BhA(t)*y(2);
end

