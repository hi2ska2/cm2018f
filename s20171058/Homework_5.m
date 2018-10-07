clear all;

k_B = 1.380648521e-23; %m^2 kg s^-2 K^-1
q = 1.602192e-19; % C
N_p(1,:) = logspace(10, 18, 161); % cm^-3
N_p(2,:) = -logspace(10, 18, 161);
n_int = 10^18; % cm^-3
T = 300; % K

V_T = k_B*T;

phi = 1;
for i =1:1000
res = N_p + n_int*exp(-q*phi/(k_B*T)) - n_int*exp(q*phi/(k_B*T));
Jaco = -n_int*q/(k_B*T)*exp(-q*phi/(k_B*T)) - n_int*q/(k_B*T)*exp(q*phi/(k_B*T));

d_phi = (-res)./Jaco;

phi = phi+d_phi;
end


%set(gca, 'XScale','log', 'YScale','log')

%% analytical solution

phi2 = (k_B*T)/q*asinh((1/2)*N_p/n_int);

%%

% difference
diff_phi = (phi2-phi)./phi2 *100;
figure(1)
plot(N_p(1,:), diff_phi(1,:), 'ko', N_p(1,:), diff_phi(2,:), 'ro')
legend('positive N^+', 'negative N^+')
set(gca,'XScale','log')
xlabel('N+ (cm^-^3)')
ylabel('phi_d_i_f_f_e_r_e_n_c_e (%)')

% total
figure(2)
plot(N_p(1,:), phi(1,:) , 'bo', N_p(1,:), phi2(1,:) , 'k-', N_p(1,:), phi(2,:), 'ro', N_p(1,:), phi2(2,:) , 'k-')
set(gca,'XScale','log')
xlabel('N+ (cm^-^3)')
ylabel('phi (V)')
legend('positive N^+ (numerical)', 'positive N^+ (analytical)', 'negative N^+ (numerical)', 'negative N^+ (analytical)')

