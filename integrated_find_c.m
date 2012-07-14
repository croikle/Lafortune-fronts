function [x, u, v] = integrated_find_c(eps,h,Z,sigma,w_star)
  c_low = 1;
  c_high = 10;
  init_pt = 1e-9;

  precision = 1e-5;

  while c_high - c_low > precision
    test_c = (c_low + c_high)/2;
    [x u v] = integrated_solve(test_c,eps,h,Z,sigma,w_star,10,0.01,[1-init_pt,1-init_pt]);
    terminus = u(end);
    if terminus < 0
      c_low = test_c;
    else
      c_high = test_c;
    end
    disp(sprintf('%f ... %f', c_low, c_high));
  end
  c = test_c
end
