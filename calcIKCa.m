% Calcium-dependent K current
% ADD COMMENTS AND DOCSTRINGS
function IKCa = calcIKCa(t, mKc, Cas, V)
    gkc_ = 8.5;
    Ek = -58;
    IKCa = zeros(1, length(t));
    for i = 1:length(t)
        mKc1 = Cas(i)/(Cas(i)+0.2);
        gKc = gkc_*mKc(i)^2*mKc1;
        IKCa(i) = gKc*(V(i)-Ek);
    end
end