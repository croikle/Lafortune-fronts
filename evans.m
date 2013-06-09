% depends on adjoint_ode.m, A_ode.m
% recursive depends: A.m, fu.m, fv.m
% also the whole integrated_* tree
% specialized to epsilon = 0 for the moment
function [result, sol] = evans(eps,h,Z,sigma,w_star)
% returns a function of lambda
  [c, front, sol] = integrated_find_c(eps,h,Z,sigma,w_star);
  % where is the front basically done
  endval = sol.x(find(sol.y(1,:) < 1e-3, 1));

  options = odeset('AbsTol',1e-9,'RelTol',1e-9);
  % maybe some precision options here

  % meet in the middle
  % later: extend before 0, shift to put middle at 0
  % t1_values = [endval, 0];
  % t2_values = [0, endval];
  % length = endval/2;

  function [value,sol1,sol2] = compute(lambda,rescale,left,right,mid,varargin)
    eigenvalue = -c/2 - sqrt(c^2 + 4*lambda)/2;

    ode1 = A_ode(c,front,lambda,h,Z,sigma,w_star);
    vector1 = rescale * [1, -(c + sqrt(c^2 + 4*lambda))/2, 0];
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
    debug = (nargin > 5);
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
