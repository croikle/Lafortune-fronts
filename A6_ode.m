function result = A6_ode(eps,c,front,eigenvalue,lambda,h,Z,sigma,w_star)
  our_A = A6(eps,c,lambda,h,Z,sigma,w_star);
  rescaled_A = @(x) our_A(x) - eigenvalue * eye(6);

  function dp = ode(t,p)
    % find the front at time t
    x = front(t);
    % use A at that x value
    dp = rescaled_A(x) * p;
  end

  result = @ode;
end
