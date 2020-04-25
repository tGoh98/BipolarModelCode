function IKCa = calcIKCa(t, mKc, Cas, V)
    % calcICa Calculates the calcium-dependent K current.
    %   Input:
    %       t - Array containing the time values for each timestep the
    %           ode used (s).
    %       mKc - Probability of m Ca-dependent K gate opening.
    %       Cas - Array containing submembranuous calcium concentration for
    %           each timestep (micromolar).
    %       V - Array containing the voltage values for each timestep (mV).
    %   Output:
    %       IKCa - Array containing the calcium-dependent K current for each 
    %           timestep (pA).
    %   Assumptions:
    %       Each of the inputs must be the same size.
    
    gkc_ = 8.5; % K Ca conductance. nS.
    Ek = -58; % K reverse potential. mV.
    
    % Calculate IKCa for each timestep.
    IKCa = zeros(1, length(t));
    for i = 1:length(t)
        mKc1 = Cas(i)/(Cas(i)+0.2);
        gKc = gkc_*mKc(i)^2*mKc1;
        IKCa(i) = gKc*(V(i)-Ek);
    end
end