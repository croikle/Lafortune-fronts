% x = [u s v]'
function result = A_zero_ode(c,lambda,h,Z,sigma,w_star)
  A = [0 1 0; lambda (-c) 0; 0 0 (lambda/c)];

  function dx = ode(t,x)
    dx = A * x;

%    dx = zeros(3,1);
%    dx(1) = x(2);
%    dx(2) = lambda * x(1) - c*x(2);
%    dx(3) = lambda/c * x(3);
  end

  result = @ode;
end
