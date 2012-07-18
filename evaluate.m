function [X, U, V] = evaluate(eps,h,Z,sigma,w_star)

  [c, sol] = integrated_find_c(eps,h,Z,sigma,w_star);
  c
  X = sol.x;
  U = sol.y(1,:)';
  V = sol.y(2,:)';
end
