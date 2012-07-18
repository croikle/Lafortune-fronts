function [X, U, V] = eval_at_c(c,eps,h,Z,sigma,w_star)
  init_pt = 1e-3;
  final_time = 10;
  sol = integrated_solve(c,eps,h,Z,sigma,w_star,final_time,[1-init_pt,1-init_pt]);

  X = [0:0.01:10];
  result = deval(X, sol);
  U = result(1,:)';
  V = result(2,:)';
end
