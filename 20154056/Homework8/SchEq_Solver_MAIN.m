%% USER INPUT
    [SysParam] = Parameter_181030(); % load the parameters (for saving the simulation environment)
    
%% LOAD PARAMETERS
Load_Parameters

%% SEMI-CLASSICAL NON-LINEAR POISSON EQUATION [initial phi(z)]
[Phi] = Nonlinear_Poisson(Nz, Ei, Intf, eps, Coef_Poi, Nacc, ni, thermal);

%% SELF-CONSISTENT LOOP
Self_consistent_Loop

%% PLOT THE RESULTS
figure(10)
plot(z*1e+9, EleDens, 'or')
xlabel('Position [nm]')
ylabel('Electron Density [#/m^3]')