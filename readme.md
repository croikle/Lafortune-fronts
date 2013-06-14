Evans function example
-
The parameters I've been using:

    eps = 0; h = 0.3; Z = 6; sigma = 0.25; w_star = 0.003;  
Set up:

    e = evans(eps, h, Z, sigma, w_star);  
Now it's a function:

    e(1+2j)  
We can plot it evaluated over a circle:

    plot(e(0.3 * exp(2*pi*j*[0:0.05:1])))  
Or a grid:

    k = e(cgrid(0:2:10, 0:2:10),10);  
    plot(k); hold on; plot(k.')  

The second parameter increases accuracy.  Values needed to reach stability, for given inputs:

    10: 10
    40: 16
    100: 22
    150: 27
    200: 30
    400: 41

    50j: 13
    100j: 17
    150j: 20

    10+10j: 10
    40+40j: 16 or 17

You can find this with e.g. `arrayfun(@(x) e(10+10j,x),7:11)` and looking for
the value where the result stops changing.


Old stuff
--

This is a small Matlab script to calculate fronts for the following system:
> U' + cU - cz = 0  
> eps*V' + cV - cz = 0  
> z' = -(1/c)(1-V)*F(hU + (1-h)V)  

where F is a mess. (2.3 in the Gordon paper, with h for (1-γ⁻¹))

This system was produced by integration of
> U" + cU' + (1-V)F(W) = 0  
> eps*V" + cV' + (1-V)F(W) = 0  

where W = hU + (1-h)V.

It shoots from near (1,1) and adjusts c to land at (0,0).

Usage:
-
    [x,u,v] = integrated_find_c(eps,h,Z,sigma,w_star);
Or with some fairly arbitrary numbers:

    [x,u,v] = integrated_find_c(0.01,0.3,6,0.25,1e-3);

Then you can do things like:

    plot(u,v);
    plot(x,[u v]);

There are a few numbers one could change in integrated_solve.m and integrated_find_c.m.
Let me know what else you'd like to do, or if there are problems.  This is far from a
masterwork yet.

Acknowledgements
--
This is part of joint work with Dr. Anna Ghazaryan, Miami University, and Dr.
Stephane Lafortune, College of Charleston.
