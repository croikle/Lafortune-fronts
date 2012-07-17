% depends on integrated_ode_eps0.m, integrated_ode.m
function sol = integrated_solve(c,eps,h,Z,sigma,w_star,L,initial)
  options = odeset('RelTol',1e-9,'AbsTol',1e-9);

  x_values = [0,L];
  if eps == 0
    ode = integrated_ode_eps0(c,h,Z,sigma,w_star);
  else
    ode = integrated_ode(c,eps,h,Z,sigma,w_star);
    initial = [initial initial(2)];
    % need initial value for z too -- take V init as hack
  end
  sol = ode45(ode, x_values, initial, options);
end
