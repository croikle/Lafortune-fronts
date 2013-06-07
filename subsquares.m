function out = subsquares(in)
    a = in(1);
    b = in(2);
    c = in(3);
    d = in(4);
    x = mean([a c]);
    y = mean([b d]);
    out = { [a b x y]; [x b c y]; [a y x d]; [x y c d] };
