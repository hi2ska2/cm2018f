function [Phi] = Nonlinear_Poisson(Nz, Ei, Intf, eps, Coef_Poi, Nacc, ni, thermal)

phi = zeros(Nz, 1);
phi(:, 1) = Ei;

%% Newton method
deviation = 1;
i = 0;
while deviation > 0.001
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
            residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + ni*exp(phi(ii, 1)/thermal))*0.5;
            Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*ni*exp(phi(ii, 1)/thermal)/thermal*0.5;
        elseif ii == Intf(2)
            residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + ni*exp(phi(ii, 1)/thermal))*0.5;
            Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*ni*exp(phi(ii, 1)/thermal)/thermal*0.5;
        else
            residue(ii, 1) = residue(ii, 1) - Coef_Poi*(Nacc + ni*exp(phi(ii, 1)/thermal));
            Jacobian(ii, ii) = Jacobian(ii, ii) - Coef_Poi*ni*exp(phi(ii, 1)/thermal)/thermal;
        end
    end
    
    update = Jacobian \ (-residue);
    phi_old = phi;
    phi = phi + update;
    deviation = max(abs(1 - phi./phi_old));
    
    i = i + 1;
    fprintf('Iteration for getting initial phi : %d\n', i)
end
Phi = phi;
end