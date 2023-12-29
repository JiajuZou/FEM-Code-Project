%Q(b) h boundary
function h_value = h_function_b(x,y)
    if (y == 0 && x >= 0 && x <= 1) || (y == 1 && x >= 0 && x <= 1)
        h_value = 0.1;
    else
        h_value = 0;
    end
end
