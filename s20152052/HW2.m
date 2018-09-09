%% homework2
clear all
clc


%% input parameter
N=5; % 50, 500; % point # 
a=5.*1e-9; % scale [m]
d=a./(N-1); % 
m=0.91;
m_0=9.109e-31; %[kg]
m_eff=m*m_0; 
h_bar=1.054e-34; % [J s]

x=zeros(N,1);  % x point
for i=2:N-1
    x(i,1)=(i-1)*d;
end
x(N,1)=a; 


%% Analytical solution
% solution: sin(kx)

for i=1:N-2
%     V_Anal(:,i)=sqrt(2./a).*sin(i*pi*x./(a));
    V_Anal(:,i)=sin(i*pi*x./(a));
    D_Anal(1,i)=(i.*pi.*h_bar)^2.*6.24e18./(2.*m_eff.*a^2); 
    k_anal(i,i)=(i.*pi./(a./1e-9))^2;
end 
    V_coeff=norm(V_Anal,2);
    norm_V=(1./V_coeff)*V_Anal;
%     norm_sort_V=sortrows(norm_V,row);
[Y,I]=sort(D_Anal,'descend');

% energy & psi reindex 
for i=1:N-2
    norm_V_reidx(:,i)=norm_V(:,I(:,i)); %% 전체 vector
    D_Anl_reidx(i,i)=Y(1,i);
end
ground_anl=min(D_Anl_reidx(D_Anl_reidx~=0));


%% Numerical solution

% operation matrix
A=zeros(N,N);
for i=2:N-1
    A(i,i)=-2;
    A(i,i-1)=1; 
    A(i,i+1)=1; 
end

% due to boundary condition, psi(0)=psi(a)=0
A1=A(2:N-1,2:N-1);
[V,D]=eig(A1); % -k^2(d^2)=D (normalized V)

E=-D.*h_bar.^2.*6.24e18./(2.*m_eff*d.^2); % energy J->E : coeff:6.24e18
ground_num=min(E(E~=0)); % energy of ground state
V1=zeros(1,N-2);
V2=[V1;V;V1]; %% 전체 vector
V_Num_prob=sqrt(conj(V2).*V2);
V_Anal_prob=sqrt(conj(norm_V_reidx).*norm_V_reidx);

E_err=(ground_anl-ground_num);


