function Il = calcIl(V)
    % calcIl Calculates the leakage current.
    %   Input:
    %       V - Array containing the voltage values for each timestep (mV).
    %   Output:
    %       Il - Array containing the leakage current for each timestep (pA).
    %   Assumptions:
    %       None.

    gl = 0.23; % Leakage conductance. nS.
    El = -21; % Leakage reverse potential. mV.

    % Compute leakage current.
    Il = gl.*(V-El);
end