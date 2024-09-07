function [polyhedrons] = voronoi_diagram_2D(Params,img)
    % Number of points to plot
    numPoints = 1000;

    % Number of superellipses to make
    n_se = length(Params);

    % Parametric angle
    t = linspace(0, 2*pi, numPoints);

    % Preallocate matrices
    Xp = zeros(n_se, numPoints);
    Yp = zeros(n_se, numPoints);

    % Calculate and transform points on all the superellipses
    for i = 1:n_se
        paramsi = Params{i};
        [x, y] = superellipse(t, paramsi.a, paramsi.b, paramsi.n);
        [Xp(i,:), Yp(i,:)] = transform(x, y, paramsi.theta, paramsi.h, paramsi.k);
    end

    % Define the bounding box
    if nargin < 2 
        bbox_min = [min(Xp(:)) - 2, min(Yp(:)) - 2]; % [x_min, y_min]
        bbox_max = [max(Xp(:)) + 2, max(Yp(:)) + 2]; % [x_max, y_max]
    elseif nargin == 2
        [height, width, ~] = size(img);
        bbox_min = [0,-height];
        bbox_max = [width,0];
    end

    % Initialize B matrices and b vectors
    B_matrices = cell(n_se, 1);
    b_vectors = cell(n_se, 1);

    % Calculate the bisectors and determine the inequalities in Bx <= b form
    for i = 1:n_se
        Bi = zeros(n_se-1,2);
        bi = zeros(n_se-1,1);
        paramsi = Params{i};
        for j = 1:n_se
            if i ~= j
                paramsj = Params{j};
                [xc1, yc1, xc2, yc2] = findClosestPoints(paramsi, paramsj);

                % Calculate the line equation Ax + By + C = 0
                A = xc2 - xc1;
                B = yc2 - yc1;
                C = 0.5 * (xc1^2 + yc1^2 - xc2^2 - yc2^2);

                % Normalize the equation
                normFactor = sqrt(A^2 + B^2);
                A = A / normFactor;
                B = B / normFactor;
                C = C / normFactor;

                % Determine the inequality for the half-space
                if A * xc1 + B * yc1 + C < 0
                    Bi = [Bi; A, B];
                    bi = [bi; -C];
                else
                    Bi = [Bi; -A, -B];
                    bi = [bi; C];
                end
            end
        end
        B_matrices{i} = Bi;
        b_vectors{i} = bi;
    end

    % Lower and upper bounds from the bounding box
    l = bbox_min';
    u = bbox_max';

    % Generate the polyhedron for each superellipse
    polyhedrons = cell(n_se, 1);
    for i = 1:n_se
        rep.B = B_matrices{i};
        rep.b = b_vectors{i};
        rep.l = l;
        rep.u = u;
        polyhedrons{i} = polyh(rep, 'h');
    end

    % Plotting the superellipses
    figure;
    hold on;
    for i = 1:n_se
         fill(Xp(i,:), Yp(i,:), 'b', 'FaceAlpha', 0.3); % Adjust FaceAlpha for transparency
    end

    % Plot the bounding box
    plot([bbox_min(1), bbox_max(1), bbox_max(1), bbox_min(1), bbox_min(1)], ...
         [bbox_min(2), bbox_min(2), bbox_max(2), bbox_max(2), bbox_min(2)], 'k--', 'LineWidth', 1.5);

    % Plot each polyhedron
    for i = 1:n_se
        plot(polyhedrons{i});
    end

    % Configure the plot
    axis equal;
    grid on;
    title('Voronoi Diagram in 2D');
    xlabel('x');
    ylabel('y');
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

function d = distanceBetweenSuperellipses(t1, t2, params1, params2)
    [x1, y1] = superellipse(t1, params1.a, params1.b, params1.n);
    [x2, y2] = superellipse(t2, params2.a, params2.b, params2.n);
    [x1, y1] = transform(x1, y1, params1.theta, params1.h, params1.k);
    [x2, y2] = transform(x2, y2, params2.theta, params2.h, params2.k);
    d = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

function [xc1, yc1, xc2, yc2] = findClosestPoints(params1, params2)
    % Objective function for optimization
    objective = @(t) distanceBetweenSuperellipses(t(1), t(2), params1, params2);

    % % Bounds for t1 and t2
    % lb = [0, 0];
    % ub = [2*pi, 2*pi];

    % Use particleswarm for optimization
    %options = optimoptions('particleswarm', 'Display', 'off');
    % Optimization options for particleswarm
    % options = optimoptions('particleswarm', ...
    %     'Display', 'off', ...           % Turn off display of iterations
    %     'UseParallel', false, ...       % Do not use parallel computation
    %     'SwarmSize', 100, ...           % Increase swarm size for more exploration
    %     'MaxIterations', 1000, ...      % Increase maximum iterations
    %     'FunctionTolerance', 1e-3);     % Tighten the tolerance for function change
    %[t_opt, ~] = particleswarm(objective, 2, lb, ub, options);

     % Initial guess for alpha1 and alpha2 based on the line joining centers
    h1_p = (params2.h - params1.h) * cos(params1.theta) + (params2.k - params1.k) * sin(params1.theta);
    k1_p = -(params2.h - params1.h) * sin(params1.theta) + (params2.k - params1.k) * cos(params1.theta);
    h2_p = (params1.h - params2.h) * cos(params2.theta) + (params1.k - params2.k) * sin(params2.theta);
    k2_p = -(params1.h - params2.h) * sin(params2.theta) + (params1.k - params2.k) * cos(params2.theta);
    tanalpha1 = abs(k1_p * params1.a / (h1_p * params1.b))^(1/params1.n) * sign(k1_p * params1.a / (h1_p * params1.b));
    tanalpha2 = abs(k2_p * params2.a / (h2_p * params2.b))^(1/params2.n) * sign(k2_p * params2.a / (h2_p * params2.b));
    alpha1 = atan(tanalpha1);
    alpha2 = atan(tanalpha2);
    % Adjust alpha1 and alpha2 if necessary
    [x1, y1] = superellipse(alpha1, params1.a, params1.b, params1.n);
    [x2, y2] = superellipse(alpha2, params2.a, params2.b, params2.n);
    if sign(x1 * h1_p + y1 * k1_p) < 0
        alpha1 = alpha1 + pi;
    end
    if sign(x2 * h2_p + y2 * k2_p) < 0
        alpha2 = alpha2 + pi;
    end

    options = optimoptions('fminunc','Display', 'off');
    t_opt = fminunc(objective, [alpha1, alpha2],options);

    % Calculate the closest points using the optimized parameters
    [x1, y1] = superellipse(t_opt(1), params1.a, params1.b, params1.n);
    [x2, y2] = superellipse(t_opt(2), params2.a, params2.b, params2.n);
    [xc1, yc1] = transform(x1, y1, params1.theta, params1.h, params1.k);
    [xc2, yc2] = transform(x2, y2, params2.theta, params2.h, params2.k);
end

