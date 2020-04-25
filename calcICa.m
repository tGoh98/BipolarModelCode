function ICa = calcICa(t, mCa, Cas, V)
    % calcICa Calculates the calcium current.
    %   Input:
    %       t - Array containing the time values for each timestep the
    %           ode used (s).
    %       mCa - Probability of m Ca gate opening.
    %       Cas - Array containing submembranuous calcium concentration for
    %           each timestep (micromolar).
    %       V - Array containing the voltage values for each timestep (mV).
    %   Output:
    %       ICa - Array containing the calcium current for each timestep
    %           (pA).
    %   Assumptions:
    %       Each of the inputs must be the same size.
    
    Cao = 2500; % Extracellular Ca concentration. micromolar.
    gca_ = 1.1; % Ca conductance. nS.
    
    % Calculate ICa for each timestep.
    hCa = @(t,V) exp(-(V-50)/11)/(exp(-(V-50)/11)+1);
    ICa = zeros(1, length(t));
    for i = 1:length(t)
       gCa = gca_*mCa(i)^4*hCa(t(i),V(i));
       ECa = 12.9*log(Cao/Cas(i));
       ICa(i) = gCa*(V(i)-ECa);
    end
end