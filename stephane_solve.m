% depends on stephane_ode_eps0.m, stephane_ode.m
function output_array = stephane_solve(c,eps,h,Z,sigma,w_star,L,dx,initial)
  options = odeset('RelTol',1e-9,'AbsTol',1e-9);

% u = [U Ux V]'
  %init_U = 0;
  %init_Ux = -0.1;
  %init_V = 0;
  %initial = [init_U init_Ux init_V];

  x_values = [dx:dx:L];
  if eps == 0
    ode = stephane_ode_eps0(c,h,Z,sigma,w_star);
    [X, Y] = ode45(ode, x_values, initial, options);
    Vx = gradient(Y(:,3),dx);
    output_array = [X Y Vx];
  else
    ode = stephane_ode(c,eps,h,Z,sigma,w_star);
    [X, Y] = ode45(ode, x_values, [initial -0.01], options);
    output_array = [X Y];
  end

  %U = Y(:,1);
  %Ux = Y(:,2);
  %V = Y(:,3);
  %Vx = gradient(Y(:,3));

  %plot(Y(:,1),Y(:,3));
  %axis([0 1 0 1]);
  %plot(Y(end,1),Y(end,3),'o');
end


