deviation1 = 1;
k = 0;
while deviation1 > 0.001
    %% Schrodinger equation, psi_zn(z), E_zn
    Phi_old = Phi;
    V = q*(EcminusEi - Phi); % Potential energy for conduction band electrons [eV]
    Nbulk = Intf(2) - Intf(1) - 1; % Number of nodes in silicon region.

    % hamiltonian
    H = zeros(Nbulk, Nbulk);
    H(1, 1) = -2; H(1, 2) = 1;
    for ii = 2 : Nbulk-1
        H(ii, ii+1) = 1;
        H(ii, ii) = -2;
        H(ii, ii-1) = 1;
    end
    H(Nbulk, Nbulk) = -2; H(Nbulk, Nbulk-1) = 1;
    for ii = 1 : Nbulk
        H(ii, ii) = H(ii, ii) - 2*m(3)*m0*(Deltaz/hbar)^2 * V(ii+Intf(1), 1);
    end

    % solution
    [eigenvectors, eigenvalues] = eig(H);
    scaledEz = diag(eigenvalues)/(-2*m(3)*m0*(Deltaz/hbar)^2); % Eigenenergy [eV]
    [sortedEz, sortedIdx] = sort(scaledEz);

    %% Calculation of the electron density n(z)
    nSubband = 10;
    EleDens = zeros(Nz, 1); % Electron density [#/m^3]

    totalEleNumber = 0;
    for n = 1 : nSubband
        Ez = sortedEz(n, 1);
        wavefunction2 = eigenvectors(:, sortedIdx(n)).^2;
        normalization = sum(wavefunction2)*Deltaz;
        wavefunction2 = wavefunction2 / normalization;
        subbandEleNumber = Coef_Sch * log(1+exp(-Ez/(k_B*T)));
        totalEleNumber = totalEleNumber + subbandEleNumber;

        EleDens(Intf(1)+1 : Intf(2)-1, 1) = EleDens(Intf(1)+1 : Intf(2)-1, 1) + 1/(L(1)*L(2))*wavefunction2*subbandEleNumber;
    end

    %% Poisson equation phi(z)
    deviation2 = 1;
    phi = Phi;
    j = 0;
    while deviation2 > 0.001
        % initializing
        residue = zeros(Nz, 1);
        Jacobian = sparse(Nz, Nz);

        % boundary condition
        residue(1, 1) = phi(1, 1) - Ei;
        Jacobian(1, 1) = 1;
        residue(Nz, 1) = phi(Nz, 1) - Ei;
        Jacobian(Nz, Nz) = 1;

        % building up matrix
        for ii = 2 : Nz-1
            if ii < Intf(1) || ii > Intf(2)
                residue(ii, 1) = eps(1)*phi(ii+1, 1) - 2*eps(1)*phi(ii, 1) + eps(1)*phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps(1); Jacobian(ii, ii) = -2*eps(1); Jacobian(ii, ii+1) = eps(1);
            elseif ii == Intf(1)
                residue(ii, 1) = eps(2)*phi(ii+1, 1) - (eps(2)+eps(1))*phi(ii, 1) + eps(1)*phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps(1); Jacobian(ii, ii) = -(eps(2)+eps(1)); Jacobian(ii, ii+1) = eps(2);
            elseif ii == Intf(2)
                residue(ii, 1) = eps(3)*phi(ii+1, 1) - (eps(3)+eps(2))*phi(ii, 1) + eps(2)*phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps(2); Jacobian(ii, ii) = -(eps(3)+eps(2)); Jacobian(ii, ii+1) = eps(3);
            else
                residue(ii, 1) = eps(2)*phi(ii+1, 1) - 2*eps(2)*phi(ii, 1) + eps(2)*phi(ii-1, 1);
                Jacobian(ii, ii-1) = eps(2); Jacobian(ii, ii) = -2*eps(2); Jacobian(ii, ii+1) = eps(2);
            end
        end

        for ii = Intf(1) : Intf(2)
            if ii == Intf(1)
                residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + EleDens(ii, 1))*0.5;
                Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*EleDens(ii, 1)/thermal*0.5;
            elseif ii == Intf(2)
                residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + EleDens(ii, 1))*0.5;
                Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*EleDens(ii, 1)/thermal*0.5;
            else
                residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + EleDens(ii, 1));
                Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*EleDens(ii, 1)/thermal;
            end
        end

        update = Jacobian \ (-residue);
        phi_old = phi;
        phi = phi + update;
        deviation2 = max(abs(1 - phi./phi_old));
        
        j = j + 1;
        fprintf('Iteration for poisson equation in self consistent loop : %d\n', j)
    end
    Phi = phi;
    deviation1 = max(abs(1 - Phi./Phi_old));
    
    k = k + 1;
    fprintf('Iteration for self consistent loop : %d\n', k)
end