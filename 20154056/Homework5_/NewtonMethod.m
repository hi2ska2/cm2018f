%% System Parameters
kB = 1.38e-23; % [J/K]
T = 300; % [K]
nint = 1e+10 * 1e+6; % [m^-3]
Nplus = 1e+15 * 1e+6; % [m^-3]

%% Numerical Solution
PhiN = [4; -4] * 1e-20; % [J]
for i = 1 : 10
    J = [-nint/(kB*T)*(exp(-PhiN(1)/(kB*T))+exp(PhiN(1)/(kB*T))), 0;
         0, -nint/(kB*T)*(exp(-PhiN(2)/(kB*T))+exp(PhiN(2)/(kB*T)))];
     
    res = [Nplus + nint*exp(-PhiN(1)/(kB*T)) - nint*exp(PhiN(1)/(kB*T));
           -Nplus + nint*exp(-PhiN(2)/(kB*T)) - nint*exp(PhiN(2)/(kB*T))];
       
    deltaPhiN = J \ (-res);
    PhiN = PhiN + deltaPhiN;
end

%% Analytical Solution
PhiA1 = kB*T*asinh(Nplus/(2*nint));
PhiA2 = kB*T*asinh(-Nplus/(2*nint));

PhiA = [PhiA1; PhiA2];

%% Display the results
PhiA
PhiN