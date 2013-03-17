% x = [u s v]'
% depends on fu.m, fv.m
function result = A(c,lambda,h,Z,sigma,w_star)
  Fu = fu(h,Z,sigma,w_star);
  Fv = fv(h,Z,sigma,w_star);
  function it=foo(x)
    it = [0, 1, 0; lambda-Fu(x), -c, -Fv(x); -Fu(x)/c, 0 (lambda-Fv(x))/c];
  end

  result = @foo;
end
