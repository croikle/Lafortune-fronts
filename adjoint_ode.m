% x = [u s v]'
function result = A_ode(c,lambda,h,Z,sigma,w_star)
  our_A = A(c,lambda,h,Z,sigma,w_star);

  function dx = ode(t,x)
    dx = our_A(x).' * x;
  end

  result = @ode;
end
