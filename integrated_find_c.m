function [c, sol] = integrated_find_c(eps,h,Z,sigma,w_star)
  c_low = 1;
  c_high = 10;
  init_pt = 1e-9;
  final_time = 10;

  precision = 1e-5;

  while c_high - c_low > precision
    test_c = (c_low + c_high)/2;
    sol = integrated_solve(test_c,eps,h,Z,sigma,w_star,final_time,[1-init_pt,1-init_pt]);
    terminus = deval(sol,final_time,1);
    if terminus < 0
      c_low = test_c;
    else
      c_high = test_c;
    end
    disp(sprintf('%f ... %f', c_low, c_high));
  end
  c = test_c;
end
