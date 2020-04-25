% Hyperpolarization activated current
function Ih = calcIh(t, O1, O2, O3, V)
    gh_ = 0.975;    
    mh = O1+O2+O3;
    Eh = -17.7;
    gh = gh_.*mh;
    Ih = zeros(1, length(t));
    for i = 1:length(t)
        Ih(i) = gh(i)*(V(i)-Eh);
    end
end