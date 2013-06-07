function out = refine(func, sq)
    found = false;
    subs = subsquares(sq);
    for i = 1 : 4
        mapped = arrayfun(func, square(subs{i}));
        if inpolygon(0,0,real(mapped),imag(mapped))
            out = subs{i};
            found = true;
            break;
        end
    end
    if ~found
        error('failed')
    end
end
