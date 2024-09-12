function [polyhedrons] = Copy_of_voronoi_diagram_2D(Params, img)
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
        bbox_min = [0, -height];
        bbox_max = [width, 0];
    end

    % Step 1: Clustering based on closest points with transitive closure (indirect overlap)
    clusters = {};
    cluster_map = zeros(1, n_se); % Cluster assignment for each superellipse
    cluster_id = 0; % Initial cluster ID

    for i = 1:n_se
        if cluster_map(i) == 0 % If not already part of a cluster
            cluster_id = cluster_id + 1; % Create a new cluster
            cluster_map(i) = cluster_id;
            clusters{cluster_id} = i;
        end

        for j = i+1:n_se
            [xi, yi, xj, yj] = findClosestPoints(Params{i}, Params{j});
            if sqrt((xi - xj)^2 + (yi - yj)^2) < 1e-3
                % If they are close enough, merge them into the same cluster
                if cluster_map(j) == 0
                    % If not yet assigned, assign it to the same cluster as i
                    cluster_map(j) = cluster_map(i);
                    clusters{cluster_map(i)} = [clusters{cluster_map(i)}, j];
                elseif cluster_map(j) ~= cluster_map(i)
                    % If they are in different clusters, merge the clusters
                    cluster_to_merge = cluster_map(j);
                    clusters{cluster_map(i)} = [clusters{cluster_map(i)}, clusters{cluster_to_merge}];
                    clusters{cluster_to_merge} = [];
                    
                    % Update the cluster map for the merged cluster
                    for k = 1:n_se
                        if cluster_map(k) == cluster_to_merge
                            cluster_map(k) = cluster_map(i);
                        end
                    end
                end
            end
        end
    end

    % Remove empty clusters
    clusters = clusters(~cellfun('isempty', clusters));

    % Step 2: Form polygons for each cluster and calculate bisectors between clusters
    num_clusters = length(clusters);
    polyhedrons = cell(num_clusters, 1);
    B_matrices = cell(num_clusters, 1);
    b_vectors = cell(num_clusters, 1);

    for i = 1:num_clusters
        cluster_i = clusters{i};
        Bi_cluster = [];
        bi_cluster = [];

        % Calculate bisectors between this cluster and other clusters
        for j = 1:num_clusters
            if i ~= j
                cluster_j = clusters{j};
                [closest_xi, closest_yi, closest_xj, closest_yj] = findClosestPointsBetweenClusters(Params, cluster_i, cluster_j);

                % Calculate bisector equation Ax + By + C = 0
                A = closest_xj - closest_xi;
                B = closest_yj - closest_yi;
                C = 0.5 * (closest_xi^2 + closest_yi^2 - closest_xj^2 - closest_yj^2);

                % Normalize
                normFactor = sqrt(A^2 + B^2);
                A = A / normFactor;
                B = B / normFactor;
                C = C / normFactor;

                % Choose center of any superellipse in cluster i to determine inequality
                params_center = Params{cluster_i(1)};
                [xc, yc] = superellipse(0, params_center.a, params_center.b, params_center.n); % Use a point on the superellipse
                [xc, yc] = transform(xc, yc, params_center.theta, params_center.h, params_center.k);

                % Determine inequality based on the side of the center
                if A * xc + B * yc + C < 0
                    Bi_cluster = [Bi_cluster; A, B];
                    bi_cluster = [bi_cluster; -C];
                else
                    Bi_cluster = [Bi_cluster; -A, -B];
                    bi_cluster = [bi_cluster; C];
                end
            end
        end

        % Bounding box inequalities for each cluster
        l = bbox_min';
        u = bbox_max';

        % Combine the bounding box inequalities
        rep.B = Bi_cluster;
        rep.b = bi_cluster;
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
    for i = 1:num_clusters
        plot(polyhedrons{i});
    end

    % Configure the plot
    axis equal;
    grid on;
    title('Voronoi Diagram with Superellipse Clusters');
    xlabel('x');
    ylabel('y');
    hold off;
end


function [closest_xi, closest_yi, closest_xj, closest_yj] = findClosestPointsBetweenClusters(Params, cluster_i, cluster_j)
    % Initialize minimum distance
    min_dist = inf;
    closest_xi = 0; closest_yi = 0; closest_xj = 0; closest_yj = 0;

    % Loop through each pair of superellipses between clusters
    for idx_i = cluster_i
        for idx_j = cluster_j
            [xi, yi, xj, yj] = findClosestPoints(Params{idx_i}, Params{idx_j});
            dist = sqrt((xi - xj)^2 + (yi - yj)^2);
            if dist < min_dist
                min_dist = dist;
                closest_xi = xi;
                closest_yi = yi;
                closest_xj = xj;
                closest_yj = yj;
            end
        end
    end
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

