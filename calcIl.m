% Leakage current
function Il = calcIl(V)
    gl = 0.23;
    El = -21;

    Il = gl.*(V-El);
end