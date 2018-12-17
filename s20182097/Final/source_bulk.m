function [source_bulk_value, source_bulk_Jaccobian] = source_bulk(phi, x_node_width, y_node_width, Dirichlet_cheak, silicon_index, Nacc_1, Nacc_2, Nacc_3, cellnumber)

q = 1.602192e-19;
ni = 1.075e16;
K_B = 1.380662e-23;
T = 300;
eps_0 = 8.8541878176e-12;


index = (1:cellnumber)';
Nacc = Nacc_1 + Nacc_2 + Nacc_3;

electron_density = x_node_width * y_node_width * q * ni * exp(q * phi /(K_B * T))/eps_0;
doping_density = x_node_width * y_node_width * q * Nacc/eps_0;
electron_density_Jacco = x_node_width * y_node_width * q * ni * (q/(K_B * T))*exp(q * phi /(K_B * T))/eps_0;

source_bulk_index = silicon_index.*double(~Dirichlet_cheak).*index;
source_bulk = [silicon_index.*double(~Dirichlet_cheak) electron_density doping_density electron_density_Jacco];

source_bulk_index( ~any(source_bulk_index, 2), : ) = [];
source_bulk( ~any(source_bulk(:,1), 2), : ) = [];

source_bulk_index_size = size(source_bulk_index);                        
source_bulk_value = sparse(source_bulk_index, ones(source_bulk_index_size(1, 1), 1), source_bulk(:, 2) + source_bulk(:, 3), cellnumber, 1);
source_bulk_Jaccobian = sparse(source_bulk_index, source_bulk_index, source_bulk(:, 4) , cellnumber, cellnumber);

end