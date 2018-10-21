h = 6.62617e-34; % Plannck constant, Js
hbar = h / (2*pi);
q = 1.602192e-19;
m0 = 9.109534e-31;
k_B = 1.308662e-23;
T = 300.0;
Lx = 100e-9; Ly = 100e-9; Lz = 5e-9; %lengths,m
mxx = 0.19; myy = 0.19; mzz = 0.91; %masses, m0
nmax = 10;
coef = 2*Lx*Ly/(2*pi)*sqrt(mxx*myy)*m0/(hbar^2)*(k_B*T);

totalNumber = 0;
Nz = 51;
dx = Lz/(Nz-1);
z = transpose([0:Nz-1])*Lz/(Nz-1);
elec = zeros(Nz,21); %electron density, /m^3
integ_n(11,1) = zeros;
Fermi(11,1) = zeros;
for ii=1:11
    for n=1:10
        Ez = (hbar^2)/(2*mzz*m0)*(pi*n/Lz)^2;
        subbandNumber = coef*log(1+exp((-Ez+((ii*2-12)*q/100))/(k_B*T)));
        totalNumber = totalNumber + subbandNumber;
        elec(:,ii) = elec(:,ii) + 2/(Lx*Ly*Lz)*(sin(n*pi*z/Lz).^2)*subbandNumber;
    end
    
    Fermi(ii,1) = (ii*2-12)/100;
    
    for jj=1:Nz
        if(jj==1 || jj==Nz)
            integ_n(ii,1) = integ_n(ii,1) + elec(jj,ii)*dx/2;
        else
            integ_n(ii,1) = integ_n(ii,1) + elec(jj,ii)*dx;
        end
    end 
    figure(1)
    plot(z/1e-9,elec(:,ii)/1e6);
    hold on;
end
xlabel('Position (nm)');
ylabel('electron density (cm^-^3)');
legend ('-0.1 eV','-0.08 eV','-0.06 eV','-0.04 eV','-0.02 eV','0 eV','0.02 eV','0.04 eV','0.06 eV','0.08 eV','0.1 eV');

figure(2)
semilogy(Fermi(:,1),integ_n(:,1)/1e4);
xlabel ('Fermi energy (eV)');
ylabel ('Integrated electron density (cm^-^2)');
