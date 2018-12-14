clear all
clc

q = 1.602192e-19; % Elementary charge, C 
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
T = 300.0; % Temperature, K 
thermal = k_B*T/q; % Thermal voltage, V 
Deltax = 1e-9; % 1 nm spacing 
N = 601; % 600-nm-long structure & 120-nm
x = Deltax*transpose([0:N-1]); % real space, m 
x_12 = 101; % At x=100 nm 40 nm
x_23 = 501; % At x=500 nm 80 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity 
Ndon = 2e21*ones(N,1); % 2e15 /cm^3 
Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3 
Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3 
ni = 1.075e16; % 1.075e10 /cm^3 
coef = Deltax*Deltax*q/eps0; 
mob_si=0.14; % Si mobility m^2/Vs
Dn_Si=thermal*mob_si;



for bias=0:10
    
    phi=zeros(N,1);
    phi(:,1)=thermal*log(Ndon(:,1)/ni);
    elec=zeros(N,1);
    elec=ni*exp(phi/thermal);
    V_applied(bias+1,1)=0.05*bias;

for newton=1:10
    
%% Poisson euqation
res=zeros(2*N,1);
Jaco=sparse(2*N,2*N);
res(1,1)=phi(1,1)-thermal*log(Ndon(1,1)/ni);
Jaco(1,1)=1.0;
  
for ii=2:N-1
    res(2*ii-1,1)=eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1))+coef*(Ndon(ii,1)-elec(ii,1));
    Jaco(2*ii-1,2*ii+1)=eps_si;
    Jaco(2*ii-1,2*ii-1)=-2*eps_si;
    Jaco(2*ii-1,2*ii-3)=eps_si;
    Jaco(2*ii-1,2*ii)=-coef;
end
    
res(2*N-1,1)=phi(N,1)-thermal*log(Ndon(N,1)/ni)-V_applied(bias+1,1);
Jaco(2*N-1,2*N-1)=1.0;


%% Continuity equation

for ii=1:N-1 % edge-wise construction    
    n_av = 0.5*(elec(ii+1,1)+elec(ii,1));   
    dphidx = (phi(ii+1,1)-phi(ii,1))/Deltax;    
    delecdx = (elec(ii+1,1)-elec(ii,1))/Deltax;    
    Jn = n_av * dphidx - thermal * delecdx;    
    res(2*ii,1) = res(2*ii,1) + Jn;    
    Jaco(2*ii,2*ii+2) = Jaco(2*ii,2*ii+2) + 0.5*dphidx - thermal / Deltax;     
    Jaco(2*ii,2*ii  ) = Jaco(2*ii,2*ii  ) + 0.5*dphidx + thermal / Deltax;      
    % electric potential에대한 의존성 포함
    Jaco(2*ii,2*ii+1)=Jaco(2*ii,2*ii+1)+n_av/Deltax;
    Jaco(2*ii,2*ii-1)=Jaco(2*ii,2*ii-1)-n_av/Deltax;
    
    res(2*ii+2,1) = res(2*ii+2,1) - Jn;    
    Jaco(2*ii+2,2*ii+2) = Jaco(2*ii+2,2*ii+2) - 0.5*dphidx + thermal / Deltax;     
    Jaco(2*ii+2,2*ii  ) = Jaco(2*ii+2,2*ii  ) - 0.5*dphidx - thermal / Deltax;      
    % electric potential에대한 의존성 포함
    Jaco(2*ii+2,2*ii+1) = Jaco(2*ii+2,2*ii+1) - n_av/Deltax;  
    Jaco(2*ii+2,2*ii-1) = Jaco(2*ii+2,2*ii-1) + n_av/Deltax;
end

% boundary condition
res(2,1) = elec(1,1) - Ndon(1,1); 
Jaco(2,:) = 0.0; 
Jaco(2,2) = 1.0;  
res(2*N,1) = elec(N,1) - Ndon(N,1); 
Jaco(2*N,:) = 0.0; 
Jaco(2*N,2*N) = 1.0;

% Scaling
Cvector=zeros(2*N,1);
Cvector(1:2:2*N-1,1)=thermal;
Cvector(2:2:2*N,1)=max(abs(Ndon));
Cmatrix=spdiags(Cvector,0,2*N,2*N); %sparse diagonal matrix로 변환
Jaco_scaled=Jaco*Cmatrix;
Rvector=1./sum(abs(Jaco_scaled),2);
Rmatrix=spdiags(Rvector,0,2*N,2*N);
Jaco_scaled=Rmatrix*Jaco_scaled;
res_scaled=Rmatrix*res;
update_scaled=Jaco_scaled\(-res_scaled);
update=Cmatrix*update_scaled;

phi=phi+update(1:2:2*N-1,1);
elec=elec+update(2:2:2*N,1);
norm(update(1:2:2*N-1,1),inf); % error monitoring

end
phi1(:,bias+1)=phi;
elec1(:,bias+1)=elec;
end

for bias=0:10
    for ii=2:N
        B=(phi1(ii,bias+1)-phi1(ii-1,bias+1))./thermal;
        J(ii-1,bias+1)=-q*Dn_Si.*(elec1(ii,bias+1)*B/(exp(B)-1)-elec1(ii-1,bias+1)*(-B)/(exp(-B)-1))/Deltax;
    end
end


figure
plot(V_applied,J(N-1,:))
axis([0 inf 0 inf])
xlabel('Applied Voltage (V)')
ylabel('Terminal Current (A/m^2)')
legend('Long structure')