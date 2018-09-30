q = 1.602e-19; % Elementary charge, C
e0 = 8.854e-12; % Vacuum permittivity, F/m
k_B = 1.381e-23; % Boltzmann constant, J/K
T = 300.0; % Temperature, K
dx = 0.1e-9; %0.1 nm spacing
N = 61; % 6 nm thick
x = dx*transpose ([0:N-1]); % real space, m
interface1 = 6; %At x=0.5 nm
interface2 = 56; % At x=5.5 nm
e_si = 11.7; e_ox = 3.9; % Relative permittivity
Nacc = 1e24; % 1e18 /cm^2
ni = 1.075e16; % 1.075e10 /cm^3

A = zeros(N,N);
A(1,1) = 1.0;
for i=2:N-1
    if     (i< interface1) A(i,i-1) = e_ox; A(i,i) = -2*e_ox;      A(i,i+1) = e_ox;
    elseif (i==interface1) A(i,i-1) = e_ox; A(i,i) = -e_ox-e_si; A(i,i+1) = e_si;
    elseif (i< interface2) A(i,i-1) = e_si; A(i,i) = -2*e_si;      A(i,i+1) = e_si;
    elseif (i==interface2) A(i,i-1) = e_si; A(i,i) = -e_ox-e_si; A(i,i+1) = e_ox;
    elseif (i> interface2) A(i,i-1) = e_ox; A(i,i) = -2*e_ox;      A(i,i+1) = e_ox;
    end
end
A(N,N) = 1.0;

b = zeros (N,11);
phi = zeros (N,11);
for j=1:11
    b(1,j) = 0.33374 - (j-1)/10.0;
    for i=interface1:interface2
	  if(i==interface1) b(i,j) = dx*dx*q*Nacc/e0*0.5;
      elseif(i==interface2) b(i,j) = dx*dx*q*Nacc/e0*0.5;
      else b(i,j)= dx*dx*q*Nacc/e0;
      end
    end
    b(N,j) = 0.33374 - (j-1)/10.0;
    phi(:,j) = A\b(:,j);
end

elec = zeros(N,11);
for j=1:11
    for i=interface1:interface2
         elec(i,j) = ni * exp(q*phi(i,j)/(k_B*T));
    end
end

c = zeros (N,11);
for j=1:11
    c(1,j) = 0.33374 - (j-1)/10.0;
    for i=interface1:interface2
        if(i==interface1) c(i,j) = dx*dx*q*(elec(i,j)+Nacc)/e0*0.5;
        elseif(i==interface2) c(i,j) = dx*dx*q*(elec(i,j)+Nacc)/e0*0.5;
        else c(i,j)= dx*dx*q*(elec(i,j)+Nacc)/e0;
        end
    end
    c(N,j) = 0.33374 - (j-1)/10.0;
    phi2(:,j)=A\c(:,j); %Updated potential
    for k=1:N
        diff(k,j)=phi(k,j)-phi2(k,j); %Potential differnce
    end
end

figure(1)
for j=1:11
    plot(x/1e-9,elec(:,j)/1e+6);
    hold on;
end
xlabel('Position (nm)');
ylabel('Electron density (cm^3)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');

figure(2)
for j=1:11
    plot(x/1e-9,phi2(:,j));
    hold on;
end
xlabel('Position (nm)');
ylabel('Updated potential (V)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');

figure(3)
for j=1:11
    plot(x/1e-9,diff(:,j));
    hold on;
end
xlabel('Position (nm)');
ylabel('Potential difference (V)');
legend('0 V','0.1 V','0.2 V','0.3 V','0.4 V','0.5 V','0.6 V','0.7 V','0.8 V','0.9 V','1.0 V','Location','bestoutside');
