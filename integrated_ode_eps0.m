% u = [U V]'
function cust_ode = integrated_ode_eps0(c,h,Z,sigma,w_star)
  function du = ode(t,u)
    du = zeros(2,1);
    du(1) = c * (u(2) - u(1));
    du(2) = -(1/c)*(1-u(2))*F(u);
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
