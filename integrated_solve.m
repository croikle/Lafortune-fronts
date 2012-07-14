% depends on integrated_ode_eps0.m, integrated_ode.m
function [X, U, V] = integrated_solve(c,eps,h,Z,sigma,w_star,L,dx,initial)
  options = odeset('RelTol',1e-9,'AbsTol',1e-9);

  x_values = [dx:dx:L];
  if eps == 0
    ode = integrated_ode_eps0(c,h,Z,sigma,w_star);
    [X, Y] = ode45(ode, x_values, initial, options);
  else
    ode = integrated_ode(c,eps,h,Z,sigma,w_star);
    [X, Y] = ode45(ode, x_values, [initial initial(2)], options);
  end
  U = Y(:,1);
  V = Y(:,2);
end
