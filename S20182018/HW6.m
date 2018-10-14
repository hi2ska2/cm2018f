q = 1.602192e-19;
eps0 = 8.854187817e-12;
k_B = 1.308662e-23;
T = 300.0;
thermal = k_B*T/q;
dx = 0.1e-9;
N = 61;
x = dx*transpose ([0:N-1]);
interface1 = 6;
interface2 = 56;
eps_si = 11.7;
eps_ox = 3.9;
Nacc = 1e24;
ni = 1.075e16;
coef = dx*dx*q/eps0;

phi = zeros(N,101);
res=zeros(N,1);
Jaco = sparse(N,N);
phi(:,1) = 0.33374;
elec = zeros(N,101);
integ_elec(101,1)=zeros;
gatev(101,1)=zeros;

for gate=1:101
    if (gate~=1)
    phi(:,gate) = phi(:,gate-1);
    end
    phi(1,gate)=0.33374-((gate-1.0)/100.0);
    phi(N,gate)=phi(1,gate);
    Jaco(1,1) = 1.0;
    Jaco(N,N) = 1.0;
    res(1,1) = phi(1,gate) - 0.33374+((gate-1.0)/100.0);
    res(N,1) = res(1,1);
    gatev(gate,1)=(gate-1.0)/100.0;
    
    for newton=1:10
        for ii=2:N-1
            if (ii< interface1 || ii>interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - 2*eps_ox*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -2*eps_ox;
                Jaco(ii,ii+1) = eps_ox;
            elseif (ii==interface1)
                res(ii,1) = eps_si*phi(ii+1,gate) - (eps_si+eps_ox)*phi(ii,gate) + eps_ox*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_ox;
                Jaco(ii,ii) = -(eps_si+eps_ox);
                Jaco(ii,ii+1) = eps_si;
            elseif (ii==interface2)
                res(ii,1) = eps_ox*phi(ii+1,gate) - (eps_ox+eps_si)*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -(eps_ox+eps_si);
                Jaco(ii,ii+1) = eps_ox;
            else
                res(ii,1) = eps_si*phi(ii+1,gate) - 2*eps_si*phi(ii,gate) + eps_si*phi(ii-1,gate);
                Jaco(ii,ii-1) = eps_si;
                Jaco(ii,ii) = -2*eps_si;
                Jaco(ii,ii+1) = eps_si;
            end
        end   
        for ii=interface1:interface2
            if (ii==interface1)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal*0.5;
            elseif(ii==interface2)
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal))*0.5;
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal*0.5;
            else
                res(ii,1) = res(ii,1) - coef*(Nacc+ni*exp(phi(ii,gate)/thermal));
                Jaco(ii,ii) = Jaco(ii,ii) - coef*ni*exp(phi(ii,gate)/thermal)/thermal;
            end
        end   
        update= Jaco \ (-res);
        phi(:,gate) = phi(:,gate) + update;
       
    end  
     for ii=interface1:interface2                    
        elec(ii,gate) = ni * exp(q*phi(ii,gate)/(k_B*T));
         if (ii== interface1 || ii==interface2)
             integ_elec(gate,1) = integ_elec(gate,1)+ dx/2 * elec(ii,gate); 
         else
             integ_elec(gate,1) = integ_elec(gate,1)+ dx * elec(ii,gate);
         end       
     end
end 

figure(1);
for i=1:101
plot (x*1e9,phi(:,i)); hold on;
end
xlabel ('Position (nm)');
ylabel ('Potential (V)');
figure(2)
for i=1:101
    plot(x*1e9,elec(:,i));hold on;
end
xlabel ('Position (nm)');
ylabel ('electron density (cm^-^3)');
figure(3)
semilogy(gatev(:,1),integ_elec(:,1)*1e-4,'o');
xlabel ('Gate voltage (V)');
ylabel ('Integrated electron density (cm^-^2)');

    
    