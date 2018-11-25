clear all
clc

q = 1.602192e-19; % Elementary charge, C 
eps0 = 8.854187817e-12; % Vacuum permittivity, F/m 
k_B = 1.380662e-23; % Boltzmann constant, J/K 
T = 300.0; % Temperature, K 
thermal = k_B*T/q; % Thermal voltage, V 
Deltax = 5e-9; % 1 nm spacing 
N = 25; % 600-nm-long structure & 120-nm
x = Deltax*transpose([0:N-1]); % real space, m 
x_12 = 9; % At x=100 nm 40 nm
x_23 = 17; % At x=500 nm 80 nm
eps_si = 11.7; eps_ox = 3.9; % Relative permittivity 
Ndon = 2e21*ones(N,1); % 2e15 /cm^3 
Ndon(1:x_12,1) = 5e23; % 5e17 /cm^3 
Ndon(x_23:N,1) = 5e23; % 5e17 /cm^3 
ni = 1.075e16; % 1.075e10 /cm^3 
coef = Deltax*Deltax*q/eps0; 



%% nonlinear poisson's equation
phi=zeros(N,1);
phi(:,1)=thermal*log(Ndon(:,1)/ni);

for newton=1:10
res=zeros(N,1);
Jaco=sparse(N,N);
res(1,1)=phi(1,1)-thermal*log(Ndon(1,1)/ni);
Jaco(1,1)=1.0;
  
for ii=2:N-1
    res(ii,1)=eps_si*(phi(ii+1,1)-2*phi(ii,1)+phi(ii-1,1));
    Jaco(ii,ii-1)=eps_si;
    Jaco(ii,ii)=-2*eps_si;
    Jaco(ii,ii+1)=eps_si;
end
    
res(N,1)=phi(N,1)-thermal*log(Ndon(N,1)/ni);
Jaco(N,N)=1.0;
    
for ii=2:N-1
    res(ii,1)=res(ii,1)-coef*(-Ndon(ii,1)+ni*exp(phi(ii,1)/thermal));
    Jaco(ii,ii)=Jaco(ii,ii)-coef*ni*exp(phi(ii,1)/thermal)/thermal;
end

update=Jaco\(-res);
phi=phi+update;
norm(update,inf);
end


%% potential 로부터 electron 구하기 위해 %%
elec=zeros(N,1);
elec=ni*exp(phi/thermal);
plot(x*1e9,elec*1e-6,'r');
hold on


%% nonlinear poisson's equation
phi=zeros(N,1);
phi(:,1)=thermal*log(Ndon(:,1)/ni);
elec=zeros(N,1);
elec=ni*exp(phi/thermal);

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
    
res(2*N-1,1)=phi(N,1)-thermal*log(Ndon(N,1)/ni);
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
Rmatrix=spdiags(Rvector,0,2*N,2*N)
Jaco_scaled=Rmatrix*Jaco_scaled;
res_scaled=Rmatrix*res;
update_scaled=Jaco_scaled\(-res_scaled);
update=Cmatrix*update_scaled;

phi=phi+update(1:2:2*N-1,1);
elec=elec+update(2:2:2*N,1);
norm(update(1:2:2*N-1,1),inf); % error monitoring
end


% 
plot(x*1e9,elec*1e-6,'o');
xlabel('Position [nm]')
ylabel('Electron density [cm^-^3]')
legend('Nonlinear Poisson equation','self0consistent')