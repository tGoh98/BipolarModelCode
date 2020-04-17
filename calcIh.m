% Hyperpolarization activated current
function [t, Ih] = calcIh(tspan, tV, voltage)
    [t,M] = ode45(@(t,M) calcProb(t,M,tV,voltage),tspan,[0;0;0;0;0]);
    mh = M(:,3)+M(:,4)+M(:,5);
    gh_ = 0.975;
    gh = gh_.*mh;
    Eh = -17.7;
    Ih = zeros(1, length(t));
    for i = 1:length(t)
        Ih(i) = gh(i)*(getV(tV,voltage,i)-Eh);
    end
end

function probabilities = calcProb(t,M,tV,voltage)
    Ah = @(t) 3/(exp((getV(tV,voltage,t)+110)/15)+1);
    Bh = @(t) 1.5/(exp(-(getV(tV,voltage,t)+115)/15)+1);
    K = @(t) -[4*Ah(t) -Bh(t) 0 0 0;
               4*Ah(t) 3*Ah(t)+Bh(t) -2*Bh(t) 0 0;
               0 -3*Ah(t) 2*Ah(t)+2*Bh(t) -3*Bh(t) 0;
               0 0 -2*Ah(t) Ah(t)+3*Bh(t) -4*Bh(t);
               0 0 0 -Ah(t) 4*Bh(t)];
    probabilities = K(t)*M;
end