function result=fv(h,Z,sigma,w_star)
  function it=fv2(x)
    if w(x) < w_star
      it = 0;
    else
      u = x(1);
      v = x(3);
      it = -(v - 1)*((sigma - 1)*(h - 1)*((h - 1)*v - h*u + h)*Z/((sigma - 1)*((h - 1)*v - h*u) + sigma)^2 - (h - 1)*Z/((sigma - 1)*((h - 1)*v - h*u) + sigma))*exp(-((h - 1)*v - h*u + h)*Z/((sigma - 1)*((h - 1)*v - h*u) + sigma)) - exp(-((h - 1)*v - h*u + h)*Z/((sigma - 1)*((h - 1)*v - h*u) + sigma));
      % this mess brought to you by computer differentiation
    end
  end
  function w=w(x)
    w = h*x(1) + (1-h)*x(3);
  end
  result = @fv2;
end
