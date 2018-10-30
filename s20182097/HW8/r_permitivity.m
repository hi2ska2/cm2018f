%% permitivity function %%
function y = r_permitivity(position, Layer_1, Layer_2, Layer_3)

    if position < Layer_1(1 , 1)
        y = Layer_1(1 , 2);
        
    elseif position < Layer_1(1 , 1) + Layer_2(1 , 1) 
        y = Layer_2(1, 2);
        
    else 
        y = Layer_3(1, 2);
        
    end
end