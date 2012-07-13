% u = [U V z]'
function cust_ode = integrated_ode(c,eps,h,Z,sigma,w_star)
  function du = ode(t,u)
    du = zeros(3,1);
    du(1) = c * (u(3) - u(1));
    du(2) = (c/eps) * (u(3) - u(2));
    du(3) = -(1/c)*(1-u(2))*F(u);
  end

  function x = F(u)
    w = h*u(1) + (1-h)*u(2);
    if w >= w_star
            x = exp(Z*(w-h)/(sigma + (1-sigma)*w));
    else
            x = 0;
    end
  end

    cust_ode = @ode;
end
