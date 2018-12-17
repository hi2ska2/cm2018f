%% Clear %%
clear all
clc
close all

%% Parameter %%
cellx = 20;
celly = 20;
cellnumber = cellx * celly;

x_length = 120e-9; %m
y_length = 6e-9; % m

x_node_width = x_length/(cellx - 1);
y_node_width = y_length/(celly - 1);

eps_0 = 8.8541878176e-12;
eps_si = 3.9;
eps_ox = 11.7; 
q = 1.602192e-19;
ni = 1.075e16;
K_B = 1.380662e-23;
T = 300;

Vg = 0;
Vd = 0;
Vs = 0;

%% Boundary & initial condition %%
index = (1:cellnumber)';

y = ceil(index/cellx);
x = index - cellx*(y - 1);

Dirichlet_cheak(index, 1) = double(  ((x_node_width * (x - 1) <= 80e-9 & x_node_width * (x - 1) >= 40e-9) & (y == 1 | y == celly))...
                                   | ((x == 1 | x == cellx) & (y_node_width * (y - 1) >= 0.5e-9 & y_node_width * (y - 1) <= 5.5e-9)));
                            
Neumann_cheak(index, 1) = double(  ((x_node_width * (x - 1) > 80e-9 | x_node_width * (x - 1) < 40e-9) & (y == 1 | y == celly))...
                                 | ((x == 1 | x == cellx) & (y_node_width * (y - 1) < 0.5e-9 | y_node_width * (y - 1) > 5.5e-9)));

silicon_index(index, 1) = double(y_node_width*(y - 1) <= 5.5e-9 & y_node_width*(y - 1) >= 0.5e-9);


Nacc_1(index, 1) = double(x_node_width*(x - 1) <= 40e-9) * 5e25; %m^-3
Nacc_2(index, 1) = double(x_node_width*(x - 1) <= 80e-9 & x_node_width*(x - 1) > 40e-9) * 2e23; %m^-3
Nacc_3(index, 1) = double(x_node_width*(x - 1) >= 80e-9) * 5e25; %m^-3


                            
%% permitivity indexing %%

permitivity_x_right(index) = eps_si * double(y_node_width * (y - 1) <= 5.5e-9 & y_node_width * (y - 1) >= 0.5e-9) + eps_ox * double(~(y_node_width * (y - 1) <= 5.5e-9 & y_node_width * (y - 1) >= 0.5e-9));
permitivity_x_left(index) = eps_si * double(y_node_width * (y - 1) <= 5.5e-9 & y_node_width * (y - 1) >= 0.5e-9) + eps_ox * double(~(y_node_width * (y - 1) <= 5.5e-9 & y_node_width * (y - 1) >= 0.5e-9));
permitivity_y_right(index) = eps_si * double(y_node_width * (y - 1 + 0.5) <= 5.5e-9 & y_node_width * (y - 1 + 0.5) >= 0.5e-9) + eps_ox * double(~(y_node_width * (y - 1 + 0.5) <= 5.5e-9 & y_node_width * (y - 1 + 0.5) >= 0.5e-9));
permitivity_y_left(index) = eps_si * double(y_node_width * (y - 1 - 0.5) <= 5.5e-9 & y_node_width * (y - 1 - 0.5) >= 0.5e-9) + eps_ox * double(~(y_node_width * (y - 1 - 0.5) <= 5.5e-9 & y_node_width * (y - 1 - 0.5) >= 0.5e-9));

permitivity = [permitivity_x_right' permitivity_x_left' permitivity_y_right' permitivity_y_left'];
                                              
%% Jaccobian %%

Jaccobian_matrix = Jaccobian(x_node_width, y_node_width, cellx, celly, Dirichlet_cheak, Neumann_cheak, permitivity);

%% V_G change 0 V to 1 V (0.1 V step)

V_G_change_phi = zeros(cellnumber, 10);
V_G_change_n = zeros(cellnumber, 10);
boundary_voltage_gate = zeros(cellnumber, 1);
boundary_voltage_drain= zeros(cellnumber, 1);
boundary_voltage_source = zeros(cellnumber, 1);

for V_G_ii = 1 : 1

    Vg = 1 + (V_G_ii - 1) * 0.1;



%% source define %%    

boundary_voltage_gate(index, 1) = (0.33374 + Vg) * double((x_node_width * (x - 1) <= 80e-9 & x_node_width * (x - 1) >= 40e-9) & (y == 1 | y == celly)); %V
boundary_voltage_drain(index, 1) = (0.33374 + Vd) * double((x == 1) & (y_node_width * (y - 1) >= 0.5e-9 & y_node_width * (y - 1) <= 5.5e-9)); %V
boundary_voltage_source(index, 1) = (0.33374 + Vs) * double((x == cellx) & (y_node_width * (y - 1) >= 0.5e-9 & y_node_width * (y - 1) <= 5.5e-9)); %V
boundary_voltage = boundary_voltage_gate + boundary_voltage_drain + boundary_voltage_source;



    phi = boundary_voltage;


Dirichlet_index = Dirichlet_cheak.*index;
Dirichlet_index( ~any(Dirichlet_index, 2), : ) = []; 
Dirichlet_index_size = size(Dirichlet_index);

boundary_voltage_merge = [Dirichlet_cheak boundary_voltage];
boundary_voltage_merge( ~any(boundary_voltage_merge(:, 1), 2), : ) = []; 

Dirichlet_value = sparse(Dirichlet_index, ones(Dirichlet_index_size(1, 1), 1), boundary_voltage_merge(:, 2), cellnumber, 1);      

[source_bulk_value, source_bulk_Jaccobian] = source_bulk(phi, x_node_width, y_node_width, Dirichlet_cheak, silicon_index, Nacc_1, Nacc_2, Nacc_3, cellnumber);
source = Dirichlet_value + source_bulk_value;





%% Non-linear poisson solver %%        
res = (Jaccobian_matrix - source_bulk_Jaccobian) * phi - source;


for Newton_ii = 1: 1000
    
    update = (Jaccobian_matrix - source_bulk_Jaccobian)\(-res);
    
    phi = phi + update;
    
    [source_bulk_value, source_bulk_Jaccobian] = source_bulk(phi, x_node_width, y_node_width, Dirichlet_cheak, silicon_index, Nacc_1, Nacc_2, Nacc_3, cellnumber);
    source = Dirichlet_value + source_bulk_value;
    res = (Jaccobian_matrix - source_bulk_Jaccobian) * phi - source;
    
    if abs(max(res)) < 5e-13
        break;
    end
    
    
end


V_G_change_phi(:, V_G_ii) = phi;
V_G_change_n(:, V_G_ii) = ni*exp(q*phi/(K_B*T)).*double(silicon_index);


end