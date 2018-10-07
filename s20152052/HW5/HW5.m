clear all
clc

q=1.602192e-19; % elementary charge, C
k_B=1.380662e-23; % Boltzmann constant
T=300; % temparature
Coeff=q/(k_B*T); % =V_T
n_int=1e10; % intrinsic carrier density of Si at 300K
N_power=(10:0.1:18)'; 
N_dop=1*10.^(N_power(:,1));% dopant density
N_pos(:,1)=1*N_dop;% dopant: donor
N_neg(:,1)=-1*N_dop;% dopant: acceptor

%% numerical solution
for i=1:length(N_power) % Donor
    phi=1; % initial value
    for newton=1:60
        Jaco=-Coeff*n_int*exp(-Coeff*phi)-Coeff*n_int*exp(Coeff*phi); % jacobian
        res=N_pos(i,1)+n_int*exp(-Coeff*phi)-n_int*exp(Coeff*phi); %residue
        update=Jaco\(-res);
        phi=phi+update;
    end
    phi_pos_num(i,1)=phi;
end

for i=1:length(N_power) % Acceptor
    phi=-1; % initial value
    for newton=1:60
        Jaco=-Coeff*n_int*exp(-Coeff*phi)-Coeff*n_int*exp(Coeff*phi); % jacobian
        res=N_neg(i,1)+n_int*exp(-Coeff*phi)-n_int*exp(Coeff*phi); %residue
        update=Jaco\(-res);
        phi=phi+update;
    end
    phi_neg_num(i,1)=phi;
end



%% analytic solution
for i=1:length(N_power)
    phi_pos_anal(i,1)=(1/Coeff)*asinh(N_pos(i,1)./(2*n_int)); % Donor
    phi_neg_anal(i,1)=(1/Coeff)*asinh(N_neg(i,1)./(2*n_int)); % Acceptor
end

err_pos=phi_pos_anal-phi_pos_num; % numerical solution- anlaytical solution (Donor)
err_neg=phi_neg_anal-phi_neg_num; % numerical solution- anlaytical solution (Acceptor)

% 
% subplot(2,1,1)
% semilogx(N_pos,err_pos,'o')
% xlabel('Donor density [cm^-^3]')
% ylabel('err [V]')
% subplot(2,1,2)
% semilogx(abs(N_neg),err_pos,'o')
% xlabel('Acceptor density [cm^-^3]')
% ylabel('err [V]')

semilogx(N_pos,phi_pos_num,'o',N_pos,phi_pos_anal,'b',abs(N_neg),phi_neg_num,'o',abs(N_neg),phi_neg_anal,'r')
xlabel('Dopant density [cm^-^3]')
ylabel('Eletrostatic potential [V]')


% semilogx(N_pos,phi_pos_num,'o')
% xlabel('Dopant density [cm^-^3]')
% ylabel('Eletrostatic potential [V]')


