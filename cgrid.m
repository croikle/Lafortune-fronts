function out = cgrid(xs, ys)
    [X, Y] = meshgrid(xs, ys);
    out = X + j * Y;
