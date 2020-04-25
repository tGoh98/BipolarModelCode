% Delayed rectifying potassium current
function IKv = calcIKv(t, mKv, hKv, V)
    gkv_ = 2.0;
    Ek = -58;
    IKv = zeros(1, length(t));
    for i = 1:length(t)
        gKv = gkv_*mKv(i)^3*hKv(i);
        IKv(i) = gKv*(V(i)-Ek);
    end
end