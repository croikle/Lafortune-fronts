% x = [u v]
% depends on fu.m, fv.m
function result = A(eps,c,lambda,h,Z,sigma,w_star)
  Fu = fu(h,Z,sigma,w_star);
  Fv = fv(h,Z,sigma,w_star);
  function it=eps_zero(x)
    it = [0, 1, 0; lambda-Fu(x), -c, -Fv(x); -Fu(x)/c, 0, (lambda-Fv(x))/c];
  end

  function it=eps_nonzero(x)
    it = [0, 1, 0, 0; lambda-Fu(x), -c, -Fv(x), 0; 0, 0, 0, 1; -Fu(x)/eps, 0, (lambda-Fv(x))/eps, -c/eps];
  end

  if eps == 0
    result = @eps_zero;
  else
    result = @eps_nonzero;
  end

end
