% g function
% whole non-zero Dirichlet BC
function g_value = g_function(x,y)
    if (x == 0 && y>= 0 && y<= 1) || (x ==1 && y>= 0 && y<= 1) || ...
        (y == 0 && x>= 0 && x<= 1) || (y ==1 && x>= 0 && x<= 1)
        g_value = 1;
    else
        g_value = 0;
    end
end
% 

