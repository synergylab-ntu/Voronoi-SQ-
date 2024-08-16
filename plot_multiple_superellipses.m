function plot_multiple_superellipses(Params,col)
    if nargin<2 
        col = 'r';
    end
    numPoints = 100;
    n_se = size(Params,1);
    % Parametric angle
    t = linspace(0, 2*pi, numPoints);

    % Preallocate matrices
    Xp = zeros(n_se, numPoints);
    Yp = zeros(n_se, numPoints);
    hold on;
    % Calculate and transform points on all the superellipses
    for i = 1:n_se
        paramsi = Params{i};
        [x, y] = superellipse(t, paramsi.a, paramsi.b, paramsi.n);
        [Xp(i,:), Yp(i,:)] = transform(x, y, paramsi.theta, paramsi.h, paramsi.k);
        fill(Xp(i,:), Yp(i,:), col);
    end
    hold off;
end

function [x, y] = superellipse(t, a, b, n)
    x = a * sign(cos(t)) .* abs(cos(t)).^n;
    y = b * sign(sin(t)) .* abs(sin(t)).^n;
end

function [x, y] = transform(x, y, theta, h, k)
    x_rot = x * cos(theta) - y * sin(theta);
    y_rot = x * sin(theta) + y * cos(theta);
    x = x_rot + h;
    y = y_rot + k;
end