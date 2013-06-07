function out=square(in)
    a = in(1);
    b = in(2);
    c = in(3);
    d = in(4);
    out = [a + b*j, a + d*j, c + d*j, c + b*j];
