function IKv = calcIKv(t, mKv, hKv, V)
    % calcIKv Calculates the delayed rectifying potassium current.
    %   Input:
    %       t - Array containing the time values for each timestep the
    %           ode used (s).
    %       mKv - Probability of m K gate opening.
    %       hKv - Probability of h K gate opening.
    %       V - Array containing the voltage values for each timestep (mV).
    %   Output:
    %       IKv - Array containing the K current for each timestep (pA).
    %   Assumptions:
    %       Each of the inputs must be the same size.

    gkv_ = 2.0; % K conductance. nS.
    Ek = -58; % K reverse potential. mV.
    
    % Calculate IKv for each timestep.
    IKv = zeros(1, length(t));
    for i = 1:length(t)
        gKv = gkv_*mKv(i)^3*hKv(i);
        IKv(i) = gKv*(V(i)-Ek);
    end
end