N0 = 1e+10;%Dopant density, cm^-3
ni = 1e+10; %intrinsic carrier density of Si, cm^-3
T = 300; %Temperature, K
k_B = 8.6173303e-5;%Boltzmann constant, eV/K
VT = k_B*T; 
phi = 10;
for ii=1:257;
    for newton=1:1000
        if (ii>=1)&&(ii<=129)
            N=-N0*10^((-ii+129)/16);
        else
            N=N0*10^((ii-129)/16);
        end
        Jaco = ni*cosh(phi/VT)/VT;
        res = -N +ni* sinh(phi/VT);
        update = Jaco \ (-res);
        phi = phi + update;
    end
    a(ii,1)=N;
    a(ii,2)=phi;
    a(ii,3)=VT*asinh(N/ni);
end
figure(1)
semilogx(a(:,1),a(:,2),'o','MarkerEdgeColor',[.4 .4 1]); hold on;
semilogx(-a(:,1),a(:,2),'o','MarkerEdgeColor',[1 .4 .4]); hold on;
semilogx(abs(a(:,1)),a(:,3),'-k','LineWidth',2);
xlabel('N^+ (cm^-^3)');
ylabel('Potential (V)');
legend('Positive N^+','Negative N^+','Exact Solution') 
figure(2)
subplot(2,1,1)
semilogx(a(:,1),(a(:,2)-a(:,3)),'.')
xlabel('Positive N^+ (cm^-^3)');
ylabel('Error (V)');
subplot(2,1,2)
semilogx(-a(:,1),(a(:,2)-a(:,3)),'.')
xlabel('Negative N^+ (cm^-^3)');
ylabel('Error (V)');
