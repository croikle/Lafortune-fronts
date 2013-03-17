% depends on adjoint_ode.m, A_ode.m
% recursive depends: A.m, fu.m, fv.m
% also the whole integrated_* tree
% specialized to epsilon = 0 for the moment
function result = evans(eps,h,Z,sigma,w_star)
% returns a function of lambda
  c = integrated_find_c(eps,h,Z,sigma,w_star);

  options = odeset('AbsTol',1e-9,'RelTol',1e-9);
  % maybe some precision options here

  endval = 5;
  t1_values = [0, -endval];
  t2_values = [0, endval];
  initial_scale = 1e-5;

  function [value,sol1,sol2] = compute(lambda)
    ode1 = A_ode(c,lambda,h,Z,sigma,w_star);
    initial1 = initial_scale*[1, -(c + sqrt(c^2 + 4*lambda))/2, 0];
    % this is the stable eigenvector of the limiting matrix at 0,0
    sol1 = ode45(ode1, t1_values, initial1, options);

    ode2 = adjoint_ode(c,lambda,h,Z,sigma,w_star);
    initial2 = initial_scale*[1, -1/2*(c + sqrt(c^2 + 4*lambda))/lambda, (c^2 + sqrt(c^2 + 4*lambda)*c)*exp(Z)/((c^2*lambda + sqrt(c^2 + 4*lambda)*c*lambda + 2*lambda^2)*exp(Z*h) + 2*lambda*exp(Z))];
    % this is the stable eigenvector of the -transposed limiting matrix at 1,1
    sol2 = ode45(ode2, t2_values, initial2, options);

    value = dot(sol1.y(:,end),sol2.y(:,end));
  end

  result = @compute;

end
