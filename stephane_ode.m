% u = [U Ux V Vx]'
function cust_ode = stephane_ode(c,eps,h,Z,sigma,w_star)
  function du = ode(t,u)
    du = zeros(4,1);
    du(1) = u(2);
    du(2) = -(c*u(2) + (1-u(3))*F(u));
    du(3) = u(4);
    du(4) = -(1/eps)*(c*u(4) + (1-u(3))*F(u));
  end

  function x = F(u)
    w = h*u(1) + (1-h)*u(3);
    if w >= w_star
            x = exp(Z*(w-h)/(sigma + (1-sigma)*w));
    else
            x = 0;
    end
  end

    cust_ode = @ode;
end
