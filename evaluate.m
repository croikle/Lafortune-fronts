function [X, U, V] = evaluate(eps,h,Z,sigma,w_star)
  X = [0:0.01:10];

  [c, sol] = integrated_find_c(eps,h,Z,sigma,w_star);
  c
  result = deval(X, sol);
  U = result(1,:)';
  V = result(2,:)';
end
