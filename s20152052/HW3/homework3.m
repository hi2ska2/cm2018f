%% homework3 : Generalized Laplacian equation

%GaAs/InAs heterostructure

clear all
clc


%% input parameter
N=101; 
a=10e-9; % scale [m]
d=a./(N-1); % 
boundary=2e-9;
e_0=8.8541878176*1e-12; % vacuum dielectric constnat (F/m)

e_1=12.4; % dielectic constant of layer1
e_2=14.6; % dielectic constant of layer2
V_l=0; % intial condition
V_r=1; % boundary condition

x=zeros(N,1);  % x point
for i=2:N-1;
    x(i,1)=(i-1)*d;
end
x(N,1)=a; 


boundary_ind=find(x==boundary)
% 빈행렬인경우: bondary or N 재설정할 것 


for i=1:N;
    if i<=boundary_ind;
        e(i,1)=e_1;
    else i>boundary_ind;
        e(i,1)=e_2;
    end
end


%% Analytical solution

syms V_anl_boundary
c=double(solve(1/(x(boundary_ind)-x(1)).*e_1.*V_anl_boundary==1/(x(N)-x(boundary_ind)).*e_2.*(V_r-V_anl_boundary))) % potential at interface
% 
syms t y
sol1=dsolve('D2y=0','y(0)=V_l','y(x(boundary_ind))=c');
sol2=dsolve('D2y=0','y(x(boundary_ind))=c','y(a)=V_r');
pretty(sol1);
pretty(sol2);
% 

V_anl=zeros(N,1); %sol1,sol2
for i=1:N;
    if i<=boundary_ind;
        V_anl(i,1)=V_l-((x(i).*(V_l-c))./x(boundary_ind)); % potential in layer1
    else i>boundary_ind;
        V_anl(i,1)=(a.*c-V_r.*x(boundary_ind))./((a-x(boundary_ind)))+(x(i).*(V_r-c)./((a-x(boundary_ind)))); % potential in layer2
    end
end

%% Numerical solution

% operation matrix
A=zeros(N,N);
A(1,1)=1;
A(N,N)=1;

for i=2:N-1;
    if i<boundary_ind;
        A(i,i)=-2.*e_1;
        A(i,i-1)=e_1; 
        A(i,i+1)=e_1;
    elseif i==boundary_ind;
        A(i,i)=-e_1-e_2;
        A(i,i-1)=e_1; 
        A(i,i+1)=e_2; 
    else i>boundary_ind;
        A(i,i)=-2.*e_2;
        A(i,i-1)=e_2; 
        A(i,i+1)=e_2;
    end
end

b=zeros(N,1); % potential
b(1,1)=V_l;
b(N,1)=V_r;


V=A\b; % potential

plot(x,V,'o',x,V_anl,'-')
xlabel('Depth (nm)')
ylabel('Potential (V)')
legend('Numerical','Analitical')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
V1_num=V(boundary_ind)-V(1);
V2_num=V(N)-V(boundary_ind);
rate_num=V1_num/V2_num;

V1_anl=V_anl(boundary_ind)-V_anl(1);
V2_anl=V_anl(N)-V_anl(boundary_ind);
rate_anl=V1_anl/V2_anl;

C1=e_1*e_0/(x(boundary_ind))
C2=e_2*e_0/(a-x(boundary_ind))
C=C1*C2/(C1+C2)

C1_num=1e-6/V1_num;
C2_num=1e-6/V2_num;
C_num=C1_num*C2_num/(C1_num+C2_num);

area=C_num/C

C1_anl=1e-6/V1_anl;
C2_anl=1e-6/V2_num;
C_anl=C1_anl*C2_anl/(C1_anl+C2_anl);
