function total = arg_principle(path)
    subtotal = 0;
    for n = 1 : (length(path) - 1)
        subtotal = subtotal + 1 - path(n)/path(n+1);
    end
    total = subtotal / (2*pi*j);
end
