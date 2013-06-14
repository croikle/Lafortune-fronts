function out = winding(point, path)
    % my_angle returns the angle between p1 and p2, as seen from point
    function a = my_angle(p1,p2)
        a1 = angle(p1-point);
        a2 = angle(p2-point);
        a = mod(a2 - a1 + pi, 2*pi) - pi;
        % returns results in [-pi,pi)
    end

    total = 0;
    for i = 1 : length(path) - 1
        total = total + my_angle(path(i), path(i+1));
    end

    out = total/(2*pi);
    
    if abs(out - round(out)) > 1e-6
        error('result not integer')
    else
        out = round(out);
    end

end
