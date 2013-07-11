% depends on adjoint_ode.m, A_ode.m
% recursive depends: A.m, fu.m, fv.m
% also the whole integrated_* tree
% specialized to epsilon = 0 for the moment
function [result, sol] = evans(eps,h,Z,sigma,w_star)
% result: a function of lambda
% sol   : the ode45 solution of the front (optional)

  [c, front, sol] = integrated_find_c(eps,h,Z,sigma,w_star);

  % Find where the front is done
  right = sol.x(find(sol.y(1,:) < 1e-2, 1))

  % Pick a point for the middle
  avg = 0.5 * (sol.y(1,:) + sol.y(2,:));
  mid = sol.x(find(avg < 0.3, 1))

  % Left is just zero
  left = 0;

  function [value,sol1,sol2] = compute_eps0(lambda, options, debug)
    eigenvalue = -c/2 - sqrt(c^2 + 4*lambda)/2;

    ode1 = A_ode(eps,c,front,lambda,h,Z,sigma,w_star);
    vector1 = [1, -(c + sqrt(c^2 + 4*lambda))/2, 0];
    % this is the stable eigenvector of the limiting matrix at 0,0
    scale1 = exp((right-mid)*eigenvalue);
    initial1 = scale1 * vector1;
    t1_values = [right, mid];
    sol1 = ode45(ode1, t1_values, initial1, options);

    ode2 = adjoint_ode(eps,c,front,lambda,h,Z,sigma,w_star);
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

  function value = compute_eps_nonzero(lambda, options, debug)
    % negative real part eigenvalues of the limit at 0,0 (+infinity)
    lambda1 = -(c + sqrt(c^2 + 4*eps*lambda))/(2*eps);
    lambda2 = -c/2 - sqrt(c^2 + 4*lambda)/2;

    % positive real part eigenvalues of the limit at 1,1 (-infinity)
    nu1 = -1/2*(c*exp(Z*h - Z) - sqrt((c^2 + 4*eps*lambda)*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z)))*exp(-Z*h + Z)/eps;
    nu2 = -c/2 + sqrt(c^2 + 4*lambda)/2;

    v1 = [0,0,1,lambda1];
    v2 = [1,lambda2,0,0];

    w1 = [1, ((sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c*eps - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c)*exp(3/2*Z*h) - (c^2*eps - c^2 - 2*eps*lambda)*exp(2*Z*h) + 2*eps*exp(Z*h + Z))/(sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps*exp(3/2*Z*h) + (2*c*eps^2 - c*eps)*exp(2*Z*h)), -1/2*(c^2*eps*exp(3*Z*h - Z) - c^2*exp(3*Z*h - Z) - ((sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*lambda*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c + 2)*eps - sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*lambda*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c)*exp(2*Z*h) + 2*(eps^2*exp(3*Z*h - Z) - eps*exp(3*Z*h - Z))*lambda)*exp(-2*Z*h)/eps^2, 1/2*(2*c*eps^2 + (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*exp(-1/2*Z*h) - 3*c)*eps - (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c^2*eps*exp(-Z*h - Z) - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*c^2*exp(-Z*h - Z) + (sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps^2*exp(-Z*h - Z) - sqrt(c^2*exp(Z*h) + 4*eps*lambda*exp(Z*h) + 4*eps*exp(Z))*eps*exp(-Z*h - Z))*lambda)*exp(3/2*Z*h) + (c^3*eps*exp(-Z*h - Z) - c^3*exp(-Z*h - Z) + 3*(c*eps^2*exp(-Z*h - Z) - c*eps*exp(-Z*h - Z))*lambda)*exp(2*Z*h))/eps^3];
% yes, it's ridiculous. This is what I got out of Sage; it may be simplifiable.
    w2 = [1,nu2,0,0];

    sigma_plus = lambda1 + lambda2;
    sigma_minus = nu1 + nu2;

    zeta_minus_prenormalized = wedge(v1, v2)
    zeta_plus = wedge(w1, w2)
    %hodge = [0 0 0 0 0 1; 0 0 0 0 -1 0; 0 0 0 1 0 0; 0 0 1 0 0 0; 0 -1 0 0 0 0; 1 0 0 0 0 0 0];
    hodge = fliplr(diag([1, -1, 1, 1, -1, 1]));

    n = dot(zeta_minus_prenormalized, hodge * zeta_plus.');
    zeta_minus = zeta_minus_prenormalized / n;

    t_values_right = [right, mid];
    t_values_left = [left, mid];

    function out = integrate(eigenvalue, eigenvector, t_values)
      len = t_values(2) - t_values(1);

      ode = A6_ode(eps,c,front,eigenvalue,lambda,h,Z,sigma,w_star);
      initial = eigenvector;
      sol = ode45(ode, t_values, initial, options);
      out = sol.y(:,end);
    end

    traceA = -c*(1 + 1/eps); % happens to be independent of \xi
    N = exp(-mid * traceA);

    left_result = integrate(sigma_minus, zeta_minus, t_values_left)
    right_result = integrate(sigma_plus, zeta_plus, t_values_right)

    value = N * dot(left_result, hodge * right_result);
  end

  function values = do_array(lambdas, varargin)
    if nargin == 1
      abstol = 6;
    else
      abstol = varargin{1};
    end

    debug = (nargin > 2);
    % have to specify abstol in order to do debug plots, sorry
    if eps == 0
      options = odeset('AbsTol',10^(-abstol),'RelTol',1e-3);
      values = arrayfun(@(x) compute_eps0(x, options, debug), lambdas);
    else
      options = odeset('AbsTol',10^(-abstol),'RelTol',10^(-abstol));
      values = arrayfun(@(x) compute_eps_nonzero(x, options, debug), lambdas);
    end
  end

  result = @do_array;
  % First argument is array of lambdas to evaluate at.
  % Second (optional) is -log10(AbsTol) (i.e. 6 -> 10^-6).
  %   For eps>0, it's applied to both absolute and relative
  % If third argument is included, plots are drawn.

end
