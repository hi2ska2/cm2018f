%% homework2 (infinite squre well)


function [ground_anl, ground_num, E_err]=E_ISW(N,param)
%% input parameter

a=param(1); % scale [m]
d=a./(N-1); 
m=param(2); 
m_0=9.109e-31; % electron mass[kg]
m_eff=m*m_0; % effective electron mass
h_bar=1.054e-34; % [J s]

A=zeros(N,N);


x=zeros(N,1);  % x point
for i=2:N-1
    x(i,1)=(i-1)*d;
end
x(N,1)=a; 


%% Analytical solution
% solution: sin(kx)

for i=1:N-2
    V_Anal(:,i)=sin(i*pi*x./(a));
    D_Anal(1,i)=(i.*pi.*h_bar)^2.*6.24e18./(2.*m_eff.*a^2); 
end 
V_coeff=norm(V_Anal,2); % normalization factor
norm_V=(1./V_coeff)*V_Anal; % normalized vector
[Y,I]=sort(D_Anal,'descend'); % sort of eivenvalue

% energy & psi reindex 
for i=1:N-2
    norm_V_reidx(:,i)=norm_V(:,I(:,i));
    D_Anl_reidx(i,i)=Y(1,i);
end
ground_anl=min(D_Anl_reidx(D_Anl_reidx~=0)); % ground energy of analitical solution


%% Numerical solution

% operation matrix
for i=2:N-1
    A(i,i)=-2;
    A(i,i-1)=1; 
    A(i,i+1)=1; 
end

% due to boundary condition, psi(0)=psi(a)=0
A1=A(2:N-1,2:N-1);
[V,D]=eig(A1); % -k^2(d^2)=D (normalized V)
E=-D.*h_bar.^2.*6.24e18./(2.*m_eff*d.^2); % J->E : coeff:6.24e18
ground_num=min(E(E~=0)); % ground energy of numerical sol.
V1=zeros(1,N-2);
V2=[V1;V;V1];
V_Num_prob=sqrt(conj(V2).*V2);
V_Anal_prob=sqrt(conj(norm_V_reidx).*norm_V_reidx);
diff=V_Num_prob(:,N-2)-V_Anal_prob(:,N-2);
plot(x,diff)

E_err=(ground_anl-ground_num);


end
