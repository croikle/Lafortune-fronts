% depends on adjoint_ode.m, A_ode.m
% recursive depends: A.m, fu.m, fv.m
% also the whole integrated_* tree
% specialized to epsilon = 0 for the moment
function [result, sol] = evans(eps,h,Z,sigma,w_star)
% result: a function of lambda
% sol   : the ode45 solution of the front (optional)

  [c, front, sol] = integrated_find_c(eps,h,Z,sigma,w_star);

  % Find where the front is done
  right = sol.x(find(sol.y(1,:) < 1e-3, 1));

  % Pick a point for the middle
  avg = 0.5 * (sol.y(1,:) + sol.y(2,:));
  mid = sol.x(find(avg < 0.3, 1));

  % Left is just zero
  left = 0;

  options = odeset('AbsTol',1e-9,'RelTol',1e-9);
  % any other options?

  function [value,sol1,sol2] = compute(lambda,varargin)
    eigenvalue = -c/2 - sqrt(c^2 + 4*lambda)/2;

    ode1 = A_ode(c,front,lambda,h,Z,sigma,w_star);
    vector1 = [1, -(c + sqrt(c^2 + 4*lambda))/2, 0];
    % this is the stable eigenvector of the limiting matrix at 0,0
    scale1 = exp((right-mid)*eigenvalue);
    initial1 = scale1 * vector1;
    t1_values = [right, mid];
    sol1 = ode45(ode1, t1_values, initial1, options);

    ode2 = adjoint_ode(c,front,lambda,h,Z,sigma,w_star);
    vector2 = [1, -1/2*(c + sqrt(c^2 + 4*lambda))/lambda, (c^2 + sqrt(c^2 + 4*lambda)*c)*exp(Z)/((c^2*lambda + sqrt(c^2 + 4*lambda)*c*lambda + 2*lambda^2)*exp(Z*h) + 2*lambda*exp(Z))];
    % this is the stable eigenvector of the -transposed limiting matrix at 1,1
    v2_normalized = vector2 / dot(vector1,vector2);
    scale2 = exp((mid-left)*eigenvalue);
    initial2 = scale2 * v2_normalized;
    t2_values = [left, mid];
    sol2 = ode45(ode2, t2_values, initial2, options);

    % for debugging. make sure these solutions don't go weird
    % pass an extra argument to show plots.
    % These solutions are actually complex, though.
    debug = (nargin > 1);
    if(debug)
      figure(1);
      plot(sol1.x,sol1.y);
      figure(2);
      plot(sol2.x,sol2.y);
    end
    value = sol1.y(:,end).' * sol2.y(:,end);
  end

  result = @compute;

end
