function Jaccobian_matrix = Jaccobian(x_node_width, y_node_width, cellx, celly, Dirichlet_cheak, Neumann_cheak, permitivity)

cellnumber = cellx*celly;
row = zeros(5*cellnumber, 1);
column = zeros(5*cellnumber, 1);
value = zeros(5*cellnumber, 1);

for ii = 1 : cellnumber
%% Boudary cheak %%

x_left = 1;
x_right= 1;
y_left = 1;
y_right = 1;

if mod(ii, cellx) == 1
    x_left = 0;
elseif mod(ii, cellx) == 0 
    x_right = 0;
end

if (mod(ii, cellx*celly) <= cellx) && (mod(ii, cellx*celly) > 0)
    y_left = 0;
elseif (mod(ii, cellx*celly) > cellx*(celly - 1)) || (mod(ii, cellx*celly) == 0)
    y_right = 0;
end


cheak = [x_left x_right y_left y_right];


%% Matrix construction %%      

        row(5*(ii - 1) + 1) = ii;
        row(5*(ii - 1) + 2) = ii;
        row(5*(ii - 1) + 3) = ii;
        row(5*(ii - 1) + 4) = ii;
        row(5*(ii - 1) + 5) = ii;


        column(5*(ii - 1) + 1) = ii + 1 * x_right;
        column(5*(ii - 1) + 2) = ii - 1 * x_left;
        column(5*(ii - 1) + 3) = ii + cellx * y_right;
        column(5*(ii - 1) + 4) = ii - cellx * y_left;
        column(5*(ii - 1) + 5) = ii;


        %Dirichlet_cheak
        if (sum(cheak) ~= 4) && (Dirichlet_cheak(ii) == 1)
            
            value(5*(ii - 1) + 1) = 0;
            value(5*(ii - 1) + 2) = 0;
            value(5*(ii - 1) + 3) = 0;
            value(5*(ii - 1) + 4) = 0;
            value(5*(ii - 1) + 5) = 1;

        %Neumann_cheak    
        elseif (sum(cheak) ~= 4) && (Neumann_cheak(ii) == 1)
            
            value(5*(ii - 1) + 1) = 1 * permitivity(ii, 1) * y_node_width/x_node_width * (0.5 +  0.5 * min(y_left, y_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * x_right;
            value(5*(ii - 1) + 2) = 1 * permitivity(ii, 2) * y_node_width/x_node_width * (0.5 +  0.5 * min(y_left, y_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * x_left;
            value(5*(ii - 1) + 3) = 1 * permitivity(ii, 3) * x_node_width/y_node_width * (0.5 +  0.5 * min(x_left, x_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * y_right;
            value(5*(ii - 1) + 4) = 1 * permitivity(ii, 4) * x_node_width/y_node_width * (0.5 +  0.5 * min(x_left, x_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * y_left;
            value(5*(ii - 1) + 5) = - (+ permitivity(ii, 1) * y_node_width/x_node_width * (0.5 +  0.5 * min(y_left, y_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * x_right... 
                                       + permitivity(ii, 2) * y_node_width/x_node_width * (0.5 +  0.5 * min(y_left, y_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * x_left...
                                       + permitivity(ii, 3) * x_node_width/y_node_width * (0.5 +  0.5 * min(x_left, x_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * y_right...
                                       + permitivity(ii, 4) * x_node_width/y_node_width * (0.5 +  0.5 * min(x_left, x_right) + 0.5 * double(abs(x_left - x_right) == 1 & abs(y_left - y_right) == 1)) * y_left);
                      
        else
            
            value(5*(ii - 1) + 1) = 1 * permitivity(ii, 1) * y_node_width/x_node_width; 
            value(5*(ii - 1) + 2) = 1 * permitivity(ii, 2) * y_node_width/x_node_width;
            value(5*(ii - 1) + 3) = 1 * permitivity(ii, 3) * x_node_width/y_node_width;
            value(5*(ii - 1) + 4) = 1 * permitivity(ii, 4) * x_node_width/y_node_width;
            value(5*(ii - 1) + 5) = -(permitivity(ii, 1) * y_node_width/x_node_width + permitivity(ii, 2) * y_node_width/x_node_width + permitivity(ii, 3) * x_node_width/y_node_width + permitivity(ii, 4) * x_node_width/y_node_width);

        end

end

Jaccobian_matrix = sparse(row, column, value, cellnumber, cellnumber);

end