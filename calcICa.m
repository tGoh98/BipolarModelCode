% Calcium current
% ADD COMMENTS AND DOCSTRINGS
function ICa = calcICa(t, mCa, Cas, V)
    Cao = 2500;
    gca_ = 1.1;
    hCa = @(t,V) exp(-(V-50)/11)/(exp(-(V-50)/11)+1);
    
    ICa = zeros(1, length(t));
    for i = 1:length(t)
       gCa = gca_*mCa(i)^4*hCa(t(i),V(i));
       ECa = 12.9*log(Cao/Cas(i));
       ICa(i) = gCa*(V(i)-ECa);
    end
end