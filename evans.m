% depends on adjoint_ode.m, A_zero_ode.m
% specialized to epsilon = 0 for the moment
function result = evans(eps,h,Z,sigma,w_star)
% returns a function of lambda
  c = integrated_find_c(eps,h,Z,sigma,w_star);

  options = odeset('AbsTol',1e-7,'RelTol',1e-7);
  % maybe some precision options here

  final = 2;
  t_values = [0,final];
  % test length, direction
  initial_scale = 1e-5;

  function [sol1,sol2] = compute(lambda)
    ode1 = A_ode(c,lambda,h,Z,sigma,w_star);
    initial1 = initial_scale*[1, -(c + sqrt(c^2 + 4*lambda))/2, 0];
    % this is the stable eigenvector of the limiting matrix at 0,0
    sol1 = ode45(ode1, -t_values, initial1, options);

    ode2 = adjoint_ode(c,lambda,h,Z,sigma,w_star);
    initial2 = [1,0,1] - initial_scale*[1, lambda/c + exp(Z*(1-h))/c, 2*lambda + (lambda^2*exp(Z*(h-1))+ exp(Z*(1-h)))/c^2];
    % this is the stable eigenvector of the limiting matrix at 1,1
    sol2 = ode45(ode2, -t_values, initial2, options);

    % evaluate them at some point
    value = {sol1 sol2};

  end

  result = @compute;

end
