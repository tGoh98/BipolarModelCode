function Ih = calcIh(t, O1, O2, O3, V)
    % calcIKv Calculates the hyperpolarization activated current.
    %   Input:
    %       t - Array containing the time values for each timestep the
    %           ode used (s).
    %       O1 - Probability of 1st open gate.
    %       O2 - Probability of 2nd open gate.
    %       O3 - Probability of 3rd open gate.
    %       V - Array containing the voltage values for each timestep (mV).
    %   Output:
    %       Ih - Array containing the hyperpolarization activated current 
    %           for each timestep (pA).
    %   Assumptions:
    %       Each of the inputs must be the same size.

    gh_ = 0.975; % Hyperpolarization conductance. nS.
    mh = O1+O2+O3; % Sum probability of open gates.
    Eh = -17.7; % Hyperpolarization reverse potential (mV).
    
    % Calculate Ih for each timestep.
    gh = gh_.*mh;
    Ih = zeros(1, length(t));
    for i = 1:length(t)
        Ih(i) = gh(i)*(V(i)-Eh);
    end
end