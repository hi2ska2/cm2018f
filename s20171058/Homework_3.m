clear all;

%%
% Permittivity
% 
% Aluminum: 1.46
% mos2: 4
% silicon dioxide: 3.9
% silicon: 11.68
%%
N = 1;

mat_th = [100 10 90 90]; %thickness
e = [1.46 4 3.9 11.68]; %permittivity

d = mat_th(1);
for n = 2:numel(mat_th)
  d = gcd(d,mat_th(n));
end

Num = mat_th.*N./d;

A=zeros(sum(Num)+1,sum(Num)+1);
A(1,1)=1;
A(sum(Num)+1,sum(Num)+1)=1;

sum_num=1;
for i = 1:numel(mat_th)
    for j = sum_num+1:sum_num+Num(i)-1
        A(j,j-1)= e(i);
        A(j,j) = -2*e(i);
        A(j,j+1)= e(i);
    end
    
    if i<numel(mat_th)
    sum_num=sum_num+Num(i);
    A(sum_num,sum_num-1) = e(i);
    A(sum_num,sum_num) = -e(i) -e(i+1);
    A(sum_num,sum_num+1) = e(i+1);
    end
end

b= zeros(sum(Num)+1,1);
b(sum(Num)+1,1) = 1;
x = A \ b;

%% analytical solution

C(1) = e(1)/mat_th(1)*10^9;
C(2) = e(2)/mat_th(2)*10^9;
C(3) = e(3)/mat_th(3)*10^9;
C(4) = e(4)/mat_th(4)*10^9;

V(1)=1/(1+C(1)/C(2)+C(1)/C(3)+C(1)/C(4));
V(2)=C(1)/C(2)*V(1);
V(3)=C(1)/C(3)*V(1);
V(4)=C(1)/C(4)*V(1);


%% plot

V_tot(1) = 0
V_tot(2) = V(1)
for i = 2:numel(V)
    V_tot(i+1) = V_tot(i)+V(i)
end

Num_x(1) = 0
Num_x(2) = Num(1)
for i = 2:numel(Num)
    Num_x(i+1) = Num_x(i)+Num(i)
end
figure(1)
xdata = 0:sum(Num)
plot(xdata*d/N,x,'bo',Num_x*d/N,V_tot,'r-');
xlabel('Thickness (nm)');
ylabel('Voltage (V)');
legend('numerical solution','analytical solution')