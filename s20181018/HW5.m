close all;
clear all;

n_int = 1e10; % Intrinsic carrier density of silicon at 300K
q = 1.602192e-19; % Elementary charge, C
k_B = 1.380662e-23; % Boltzmann constant
T = 300; % Temparature
VT = (k_B*T)/q;
dop_sweep = transpose(10:0.1:18);
dop = 1*10.^(dop_sweep(:,1)); % Dopant density

% Numerical Solution %
N(:,1) = 1*dop; % Donor
for i = 1:length(dop_sweep)
    phi = 1;
    for newton = 1:1000
        Jaco = n_int*(-1/VT)*exp(-phi/VT)-n_int*(1/VT)*exp(phi/VT); % Jacobian
        res = N(i,1)+n_int*exp(-phi/VT)-n_int*exp(phi/VT); % Residue
        update = Jaco\(-res);
        phi = phi+update;
    end
    phi_num(i,1) = phi;
end

N(:,2) = -1*dop; % Acceptor
for i = 1:length(dop_sweep)
    phi = -1;
    for newton=1:1000
        Jaco = n_int*(-1/VT)*exp(-phi/VT)-n_int*(1/VT)*exp(phi/VT); % Jacobian
        res = N(i,2)+n_int*exp(-phi/VT)-n_int*exp(phi/VT); % Residue
        update = Jaco\(-res);
        phi = phi+update;
    end
    phi_num(i,2) = phi;
end

% Analytic Solution %
for i = 1:length(dop_sweep)
    phi_anal(i,1) = VT*asinh(N(i,1)./(2*n_int));
    phi_anal(i,2) = VT*asinh(N(i,2)./(2*n_int));
end

diff(:,1) = phi_anal(:,1)-phi_num(:,1);
diff(:,2) = phi_anal(:,2)-phi_num(:,2);

figure(1);
semilogx(N(:,1),phi_num(:,1),'bo',N(:,1),phi_anal(:,1),'b-',abs(N(:,2)),phi_num(:,2),'ro',abs(N(:,2)),phi_anal(:,2),'r-');
xlabel('Doping concentration (cm^-^3)');
ylabel('Eletrostatic potential (V)');
legend('positive N^+ (numerical)', 'positive N^+ (analytical)', 'negative N^+ (numerical)', 'negative N^+ (analytical)','Location','Best');
figure(2);
subplot(2,1,1);
semilogx(N(:,1),diff(:,1),'bo');
xlabel('Donor concentration (cm^-^3)');
ylabel('Error (V)');
subplot(2,1,2);
semilogx(abs(N(:,2)),diff(:,2),'ro');
xlabel('Acceptor concentration (cm^-^3)');
ylabel('Error (V)');
