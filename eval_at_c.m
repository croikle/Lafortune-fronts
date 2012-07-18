function [X, U, V] = eval_at_c(c,eps,h,Z,sigma,w_star)
  init_pt = 1e-3;
  sol = integrated_solve(c,eps,h,Z,sigma,w_star,[1-init_pt,1-init_pt]);

  X = sol.x;
  U = sol.y(1,:)';
  V = sol.y(2,:)';
end
