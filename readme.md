This is a small Matlab script to calculate fronts for the following system:
> U' + cU - cz = 0  
> eps*V' + cV - cz = 0  
> z' = -(1/c)(1-V)*F(hU + (1-h)V)  
where F is a mess. (2.3 in the Gordon paper, with h for (1-gamma^(-1)))

This system was produced by integration of
> U" + cU' + (1-V)F(W) = 0  
> eps*V" + cV' + (1-V)F(W) = 0  
where W = hU + (1-h)V

It shoots from near (1,1) and adjusts c to land at (0,0).

Usage:
------
    [x,u,v] = integrated_find_c(eps,h,Z,sigma,w_star);
Or with some fairly arbitrary numbers:

    [x,u,v] = integrated_find_c(0.01,0.3,6,0.25,1e-3);

Then you can do things like:

    plot(u,v);
    plot(x,[u v]);

There are a few numbers one could change in integrated_solve.m and integrated_find_c.m.
Let me know what else you'd like to do, or if there are problems.  This is far from a
masterwork yet.
