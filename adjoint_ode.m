% x = [u v]
% p = [p s q]
% depends on A.m
% in turn: fu.m, fv.m
function result = A_ode(eps,c,front,lambda,h,Z,sigma,w_star)
  our_A = A(eps,c,lambda,h,Z,sigma,w_star);

  function dp = ode(t,p)
    % find the front at time t
    x = front(t);
    % use A at that x value
    dp = -our_A(x).' * p;
  end

  result = @ode;
end
