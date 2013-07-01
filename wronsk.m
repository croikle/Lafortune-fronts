function [out, xs] = wronsk(sols)
% sols is array of 4 ode solution objects
% all over same x values
    start_x = sols(1).x(1);
    end_x = sols(1).x(end);
    ends = arrayfun(@(s) s.x(end), sols);
    if any(ends ~= end_x)
        error 'solutions don''t end together';
    end

    % 100 points, arbitrarily
    xs = linspace(start_x, end_x);

    inp = arrayfun(@(s) deval(s, xs), sols, 'UniformOutput', false);
    % inp is 4x cell array of 4 x n arrays
    n = 100;
    for j = 1:n
        mid = zeros(4);
        for i = 1:4
            it = inp{i};
            mid(:,i) = it(:,j);
            % probably could do something else here
        end
        out(j) = det(mid);
    end
end
