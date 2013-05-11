function [c, front, sol] = integrated_find_c(eps,h,Z,sigma,w_star)
  c_low = 0;
  c_high = 10;
  initial_scale = 1e-6;

  function v = evector(c)
    if eps == 0
      v = [1, 1 + exp((1-h)*Z)/c^2];
    else
      v = [1, -1/2*(c^2*eps*exp(3*Z*h - Z) - c^2*exp(3*Z*h - Z) - ((sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c + 2)*eps - sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z))*c)*exp(2*Z*h))*exp(-2*Z*h)/eps^2, 1/2*((2*c*eps - c)*exp(Z*h) + sqrt(c^2*exp(Z*h) + 4*eps*exp(Z))*exp(1/2*Z*h))*exp(-Z*h)/(c*eps)];
      % this may be possible to simplify
    end
  end

  function v = evalue(c)
    if eps == 0
      v = exp((1-h)*Z)/c;
    else
      v = -1/2*(c*exp(Z*h - Z) - sqrt(c^2*exp(2*Z*h - 2*Z) + 4*eps*exp(Z*h - Z)))*exp(-Z*h + Z)/eps
      % this may be possible to simplify
    end
  end

  precision = 1e-6;

  if eps == 0
    fixpt = [1,1];
  else
    fixpt = [1,1,1];
  end

  while c_high - c_low > precision
    test_c = (c_low + c_high)/2;
    initial = fixpt - initial_scale * evector(test_c);
    sol = integrated_solve(test_c,eps,h,Z,sigma,w_star,initial);
    terminus = sol.y(1,end);
    if terminus < 0
      c_low = test_c;
    else
      c_high = test_c;
    end
    disp(sprintf('%f ... %f', c_low, c_high));
  end

  c = test_c;

  function y = eval_front(xi)
  if xi >= 0
      y = deval(sol, xi)'; % if xi > endpoint this fails. handle somehow?
    else
      y = fixpt - exp(xi*evalue(c)) * initial_scale * evector(c);
    end
  end
  front = @eval_front;
end
