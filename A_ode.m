% x = [u v]
% p = [p s q]
% depends on A.m
% in turn: fu.m, fv.m
function result = A_ode(c,front,lambda,h,Z,sigma,w_star)
  our_A = A(c,lambda,h,Z,sigma,w_star);
  % find "midpoint" of front
  % maybe where average of u,v = 1/2?

  function dp = ode(t,p)
    % find the front at time t
    x = deval(front,t);
    % use A at that x value
    dp = our_A(x) * p;
  end

  result = @ode;
end
