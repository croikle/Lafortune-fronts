% depends on integrated_ode_eps0.m, integrated_ode.m
function sol = integrated_solve(c,eps,h,Z,sigma,w_star,initial)
  options = odeset('RelTol',1e-9,'AbsTol',1e-9);
  final = 7;
  x_values = [0,final];

  if eps == 0
    ode = integrated_ode_eps0(c,h,Z,sigma,w_star);
  else
    ode = integrated_ode(c,eps,h,Z,sigma,w_star);
    % caller should pass in appropriate-dimension initial condition
  end
  sol = ode45(ode, x_values, initial, options);
  while abs(get_last_du(sol)) > 1e-7
    if final > 30
      disp('giving up');
      break;
      % if c is too high decline is very slow
      % assume we'd converge by now if that's not the case
    end
    final = final * 2;
    disp(['extending to ', num2str(final)]);
    sol = odextend(sol, [], final);
  end
end

function du = get_last_du(sol)
  dx = sol.x(end) - sol.x(end-1);
  dy = sol.y(1,end) - sol.y(1,end-1);
  du = dy/dx;
end
