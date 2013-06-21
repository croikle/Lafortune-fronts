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

  function [value,sol1,sol2] = compute_eps0(lambda, options, debug)
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
    if debug
      figure(1);
      plot(sol1.x,sol1.y);
      figure(2);
      plot(sol2.x,sol2.y);
    end
    value = sol1.y(:,end).' * sol2.y(:,end);
  end

  function [value,sol1,sol2,sol3,sol4] = compute_eps_nonzero(lambda, options, debug)
    % negative real part eigenvalues of the limit at 0,0 (+infinity)
    eigenvalue1 = -(c + sqrt(c^2 + 4*eps*lambda))/(2*eps);
    eigenvalue2 = -c/2 - sqrt(c^2 + 4*lambda)/2;
    % positive real part eigenvalues of the limit at 1,1 (-infinity)
    eigenvalue3 = -1/2*(c*exp(Z*h - Z) - sqrt((c^2 + 4*eps*lambda)*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z)))*exp(-Z*h + Z)/eps;
    eigenvalue4 = -c/2 + sqrt(c^2 + 4*lambda)/2;

    eigenvector1 = [0,0,1,eigenvalue1];
    eigenvector2 = [1,eigenvalue2,0,0];
    eigenvector3 = [1, ((sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c*eps - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c)*exp(3/2*Z*h) - (c^2*eps - c^2 - 2*eps*lambda)*exp(2*Z*h) + 2*eps*exp(Z*h + Z))/(sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps*exp(3/2*Z*h) + (2*c*eps^2 - c*eps)*exp(2*Z*h)), -1/2*(c^2*eps*exp(3*Z*h - Z) - c^2*exp(3*Z*h - Z) - ((sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*lambda*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c + 2)*eps - sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*lambda*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c)*exp(2*Z*h) + 2*(eps^2*exp(3*Z*h - Z) - eps*exp(3*Z*h - Z))*lambda)*exp(-2*Z*h)/eps^2, 1/2*(2*c*eps^2 + (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*exp(-1/2*Z*h) - 3*c)*eps - (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c^2*eps*exp(-Z*h - Z) - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c^2*exp(-Z*h - Z) + (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps^2*exp(-Z*h - Z) - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps*exp(-Z*h - Z))*lambda)*exp(3/2*Z*h) + (c^3*eps*exp(-Z*h - Z) - c^3*exp(-Z*h - Z) + 3*(c*eps^2*exp(-Z*h - Z) - c*eps*exp(-Z*h - Z))*lambda)*exp(2*Z*h))/eps^3];
% yes, it's ridiculous. This is what I got out of Sage; it may be simplifiable.
    eigenvector4 = [1,eigenvalue4,0,0];

    ode = A_ode(c,front,lambda,h,Z,sigma,w_star);

    scale1 = exp((right-mid)*eigenvalue1);
    initial1 = scale1 * eigenvector1;
    t1_values = [right, mid];
    sol1 = ode45(ode, t1_values, initial1, options);

    scale2 = exp((right-mid)*eigenvalue2);
    initial2 = scale2 * eigenvector2;
    t2_values = [right, mid];
    sol2 = ode45(ode, t2_values, initial2, options);
% any normalizing the scale of the eigenvectors?
% last time we ensured v1 . v2 = 1

    scale3 = exp((mid-left)*eigenvalue3);
    initial3 = scale3 * eigenvector3;
    t3_values = [left, mid];
    sol3 = ode45(ode3, t3_values, initial3, options);

    scale4 = exp((mid-left)*eigenvalue4);
    initial4 = scale4 * eigenvector4;
    t4_values = [left, mid];
    sol4 = ode45(ode4, t4_values, initial4, options);

    % for debugging. make sure these solutions don't go weird
    % pass an extra argument to show plots.
    % These solutions are actually complex, though.
    if debug
      figure(1);
      plot(sol1.x,sol1.y);
      figure(2);
      plot(sol2.x,sol2.y);
      figure(3);
      plot(sol3.x,sol3.y);
      figure(4);
      plot(sol4.x,sol4.y);
    end
    value = det([sol1.y(:,end),sol2.y(:,end),sol3.y(:,end),sol4.y(:,end)]);
  end

  function values = do_array(lambdas, varargin)
    if nargin == 1
      abstol = 6;
    else
      abstol = varargin{1};
    end
    options = odeset('AbsTol',10^(-abstol),'RelTol',1e-3);

    debug = (nargin > 2);
    % have to specify abstol in order to do debug plots, sorry
    if eps == 0
      values = arrayfun(@(x) compute_eps0(x, options, debug), lambdas);
    else
      values = arrayfun(@(x) compute_eps_nonzero(x, options, debug), lambdas);
    end
  end

  result = @do_array;
  % First argument is array of lambdas to evaluate at.
  % Second (optional) is -log10(AbsTol) (i.e. 6 -> 10^-6).
  % If third argument is included, plots are drawn.

end
