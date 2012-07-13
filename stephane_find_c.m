function c = stephane_find_c(h,Z,sigma,w_star)
  c_low = 1;
  c_high = 2;
  init_pt = 1e-4;

  precision = 1e-7;

  while c_high - c_low > precision
    test_c = (c_low + c_high)/2;
    solution = stephane_solve(test_c,h,Z,sigma,w_star,-10,-0.01,[init_pt/test_c,-init_pt,0]);
    terminus = solution(end,2);
    if terminus < 1
      c_low = test_c
    else
      c_high = test_c
    end
  end
  c = c_low;
end
