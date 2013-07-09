function result = A6_ode(eps,c,front,lambda,h,Z,sigma,w_star)
  our_A = A6(eps,c,lambda,h,Z,sigma,w_star);

  function dp = ode(t,p)
    % find the front at time t
    x = front(t);
    % use A at that x value
    dp = our_A(x) * p;
  end

  result = @ode;
end
