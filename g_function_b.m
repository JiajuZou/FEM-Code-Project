%Q(b) g boundary
function g_value = g_function_b(x,y)
    if (x == 0 && y >= 0 && y <= 1) || (x == 1 && y >= 0 && y <= 1)
        g_value = 0.1;
    else
        g_value = 0;
    end
end
