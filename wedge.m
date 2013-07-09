% wedge product C4 x C4 -> 6-d space
function out = wedge(a, b)
    out(1) = a(1) * b(2) - a(2) * b(1);
    out(2) = a(1) * b(3) - a(3) * b(1);
    out(3) = a(1) * b(4) - a(4) * b(1);
    out(4) = a(2) * b(3) - a(3) * b(2);
    out(5) = a(2) * b(4) - a(4) * b(2);
    out(6) = a(3) * b(4) - a(4) * b(3);
end
