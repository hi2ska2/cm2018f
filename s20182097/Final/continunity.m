function continuity_matrix = continunity(phi, elec, x_node_width, y_node_width, cellx, celly )

cellnumber = cellx*celly;
row = zeros(6*cellnumber, 1);
column = zeros(6*cellnumber, 1);
value = zeros(6*cellnumber, 1);

row_2nd = zeros(6*cellnumber, 1);
column_2nd = zeros(6*cellnumber, 1);
value_2nd = zeros(6*cellnumber, 1);

row_3rd = zeros(cellnumber, 1);
column_3rd = zeros(cellnumber, 1);
value_3rd = zeros(cellnumber, 1);

res = zeros(2*cellnumber, 1);

q = 1.602192e-19;
K_B = 1.380662e-23;
T = 300;
thermal = K_B*T/q;



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


%% Value %% 

if sum(cheak) == 4

    n_av_x = 0.5 * (elec(ii + 1, 1) + elec(ii, 1));
    n_av_y = 0.5 * (elec(ii + 1 * cellx, 1) + elec(ii, 1)) ;

    dphidx = (phi(ii + 1, 1) - phi(ii, 1))/x_node_width;
    dphidy = (phi(ii + 1, 1) - phi(ii, 1))/y_node_width;

    delecdx = (elec(ii + 1, 1) - elec(ii, 1))/x_node_width;
    delecdy = (elec(ii + 1, 1) - elec(ii, 1))/y_node_width;

    Jn = n_av_x * dphidx * y_node_width + n_av_y * dphidy * x_node_width - thermal * delecdx * y_node_width - thermal * delecdy * x_node_width;

else
    
    n_av_x = 0;
    n_av_y = 0;

    dphidx = 0;
    dphidy = 0;

    Jn = 0;
    
end

%% Matrix construction %%      

        row(6*(ii - 1) + 1) = 2 * ii;
        row(6*(ii - 1) + 2) = 2 * ii;
        row(6*(ii - 1) + 3) = 2 * ii;
        row(6*(ii - 1) + 4) = 2 * ii;
        row(6*(ii - 1) + 5) = 2 * ii;
        row(6*(ii - 1) + 6) = 2 * ii;
        
        column(6*(ii - 1) + 1) = 2 * (ii + 1 * x_right);
        column(6*(ii - 1) + 2) = 2 * (ii + 1 * x_right) - 1;
        column(6*(ii - 1) + 3) = 2 * (ii + cellx * y_right);
        column(6*(ii - 1) + 4) = 2 * (ii + cellx * y_right) - 1;
        column(6*(ii - 1) + 5) = 2 * ii;
        column(6*(ii - 1) + 6) = 2 * ii - 1;
       
        row_2nd(6*(ii - 1) + 1) = 2 * (ii + min(x_right, y_right));
        row_2nd(6*(ii - 1) + 2) = 2 * (ii + min(x_right, y_right));
        row_2nd(6*(ii - 1) + 3) = 2 * (ii + min(x_right, y_right));
        row_2nd(6*(ii - 1) + 4) = 2 * (ii + min(x_right, y_right));
        row_2nd(6*(ii - 1) + 5) = 2 * (ii + min(x_right, y_right));
        row_2nd(6*(ii - 1) + 6) = 2 * (ii + min(x_right, y_right));

        column_2nd(6*(ii - 1) + 1) = 2 * (ii + 1 * x_right);
        column_2nd(6*(ii - 1) + 2) = 2 * (ii + 1 * x_right) - 1;
        column_2nd(6*(ii - 1) + 3) = 2 * (ii + cellx * y_right);
        column_2nd(6*(ii - 1) + 4) = 2 * (ii + cellx * y_right) - 1;
        column_2nd(6*(ii - 1) + 5) = 2 * ii;
        column_2nd(6*(ii - 1) + 6) = 2 * ii - 1;

        row_3rd(ii, 1) = 2 * ii - 1;
        column_3rd(ii, 1) = 2 * ii; 
        value_3rd(ii, 1) = x_node_width*y_node_width*q/



               %Dirichlet_cheak
     
        if sum(cheak) ~= 4
            
            value(6*(ii - 1) + 1) = 0; 
            value(6*(ii - 1) + 2) = 0;       
            value(6*(ii - 1) + 3) = 0; 
            value(6*(ii - 1) + 4) = 0; 
            value(6*(ii - 1) + 5) = 1; 
            value(6*(ii - 1) + 6) = 0;  

        else
        
        value(6*(ii - 1) + 1) = + 0.5*dphidx * y_node_width - thermal * y_node_width/x_node_width * x_right; 
        value(6*(ii - 1) + 2) = + n_av_x * y_node_width/x_node_width * x_right;       
        value(6*(ii - 1) + 3) = + 0.5*dphidy * x_node_width - thermal * x_node_width/y_node_width * y_right; 
        value(6*(ii - 1) + 4) = + n_av_y * x_node_width/y_node_width * y_right; 
        value(6*(ii - 1) + 5) = + 0.5*dphidx * y_node_width + thermal * y_node_width/x_node_width * x_right + 0.5*dphidx * x_node_width + thermal * x_node_width/y_node_width * y_right; 
        value(6*(ii - 1) + 6) = + n_av_x * y_node_width/x_node_width * x_right + n_av_y * x_node_width/y_node_width * y_right;  
 
        value_2nd(6*(ii - 1) + 1) =  - 0.5*dphidx * y_node_width - thermal * y_node_width/x_node_width * x_right;     
        value_2nd(6*(ii - 1) + 2) = - n_av_x * y_node_width/x_node_width * x_right;       
        value_2nd(6*(ii - 1) + 3) = - 0.5*dphidy * x_node_width - thermal * x_node_width/y_node_width * y_right; 
        value_2nd(6*(ii - 1) + 4) = - n_av_y * x_node_width/y_node_width * y_right; 
        value_2nd(6*(ii - 1) + 5) = - 0.5*dphidx * y_node_width + thermal * y_node_width/x_node_width * x_right + 0.5*dphidx * x_node_width + thermal * x_node_width/y_node_width * y_right; 
        value_2nd(6*(ii - 1) + 6) = - n_av_x * y_node_width/x_node_width * x_right + n_av_y * x_node_width/y_node_width * y_right;  
        
        
        end

 %% residual construction %%             

        res(2*ii, 1) = res(2*ii, 1) + Jn;
        
        if x_right == 1 || y_right == 1 
            res(2*ii + 2, 1) = res(2*ii + 2 , 1) - Jn;
        end
end

continuity_matrix_1st = sparse(row, column, value, 2*cellnumber, 2*cellnumber);
continuity_matrix_2nd = sparse(row_2nd, column_2nd, value_2nd, 2*cellnumber, 2*cellnumber);
continuity_matrix = continuity_matrix_1st + continuity_matrix_2nd;


end