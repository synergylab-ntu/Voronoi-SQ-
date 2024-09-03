function [B_all, b_all, polyhedrons, bbox_min, bbox_max] = calculate_bisecting_planes_copy(params)
    num_superquadrics = numel(params);
    B_all = cell(num_superquadrics, 1);
    b_all = cell(num_superquadrics, 1);
    centers = zeros(num_superquadrics, 3);

    % Calculate the centers of all superquadrics
    for i = 1:num_superquadrics
        centers(i, :) = params{i}(9:11);
    end

    % Define bounding box limits
    bbox_min = min(centers) - 3;
    bbox_max = max(centers) + 3;
    l = bbox_min';
    u = bbox_max';

    % Initialize matrices to store bisecting plane results
    B_matrix = cell(num_superquadrics, num_superquadrics);
    b_matrix = cell(num_superquadrics, num_superquadrics);

    % Calculate B and b matrices for each superquadric
    for i = 1:num_superquadrics
        B_i = zeros(num_superquadrics-1, 3);
        b_i = zeros(num_superquadrics-1, 1);
        cnt = 1;
        for j = 1:num_superquadrics
            if i ~= j
                if isempty(B_matrix{i, j})
                    [u1_opt, v1_opt, u2_opt, v2_opt] = find_closest_points(params{i}, params{j});
                    [X1_opt, Y1_opt, Z1_opt] = superquadric_point(params{i}, u1_opt, v1_opt);
                    [X2_opt, Y2_opt, Z2_opt] = superquadric_point(params{j}, u2_opt, v2_opt);

                    midpoint = [(X1_opt + X2_opt) / 2, (Y1_opt + Y2_opt) / 2, (Z1_opt + Z2_opt) / 2];
                    normal = [X2_opt - X1_opt, Y2_opt - Y1_opt, Z2_opt - Z1_opt];
                    normal = normal / norm(normal); % Normalize the normal vector

                    B_matrix{i, j} = normal;
                    b_matrix{i, j} = dot(normal, midpoint);
                    B_matrix{j, i} = -normal;
                    b_matrix{j, i} = -dot(normal, midpoint);
                end
                
                B_i(cnt, :) = B_matrix{i, j};
                b_i(cnt) = b_matrix{i, j};
                cnt = cnt + 1;
            end
        end
        B_all{i} = B_i;
        b_all{i} = b_i;
    end

    % Generate the polyhedron for each superquadric
    polyhedrons = cell(num_superquadrics, 1);
    for i = 1:num_superquadrics
        rep.B = B_all{i};
        rep.b = b_all{i};
        rep.l = l;
        rep.u = u;
        polyhedrons{i} = polyh(rep, 'h');
    end

    % Plot each polyhedron
    figure;
    hold on;
    for i = 1:num_superquadrics
        plot(polyhedrons{i});
    end

    % Plot the bounding box
    plot_bounding_box(bbox_min, bbox_max);
    % Plot superquadrics
    plot_multiple_superquadrics(params,0)
    view(3);
    axis equal;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Polyhedrons and Bounding Box');
    grid on;
    hold off;
end

function plot_bounding_box(bbox_min, bbox_max)
    % Extract coordinates of the bounding box
    x_min = bbox_min(1); x_max = bbox_max(1);
    y_min = bbox_min(2); y_max = bbox_max(2);
    z_min = bbox_min(3); z_max = bbox_max(3);

    % Define vertices of the bounding box
    vertices = [
        x_min, y_min, z_min;
        x_max, y_min, z_min;
        x_max, y_max, z_min;
        x_min, y_max, z_min;
        x_min, y_min, z_max;
        x_max, y_min, z_max;
        x_max, y_max, z_max;
        x_min, y_max, z_max;
    ];

    % Define edges of the bounding box
    edges = [
        1, 2; 2, 3; 3, 4; 4, 1;
        5, 6; 6, 7; 7, 8; 8, 5;
        1, 5; 2, 6; 3, 7; 4, 8;
    ];

    % Plot the bounding box
    for i = 1:size(edges, 1)
        plot3(vertices(edges(i, :), 1), vertices(edges(i, :), 2), vertices(edges(i, :), 3), 'k-', 'LineWidth', 1.5);
    end
end

function [u1_opt, v1_opt, u2_opt, v2_opt, min_distance] = find_closest_points(params1, params2)
    % Define the distance function to be minimized
    distance_func = @(uv) compute_surface_distance(params1, params2, uv(1), uv(2), uv(3), uv(4));

    % Initial guess for u and v within bounds
    initial_guess = [0, 0, 0, 0];

    % Lower and upper bounds for u and v
    lb = [-pi/2, -pi, -pi/2, -pi];
    ub = [pi/2, pi, pi/2, pi];

    % Optimization options for particleswarm
    options = optimoptions('particleswarm', ...
        'Display', 'off', ...           % Turn off display of iterations
        'UseParallel', false, ...       % Do not use parallel computation
        'SwarmSize', 200, ...           % Increase swarm size for more exploration
        'MaxIterations', 2000, ...      % Increase maximum iterations
        'FunctionTolerance', 1e-6);     % Tighten the tolerance for function change

    % Perform optimization to find optimal u and v
    [uv_opt, min_distance] = particleswarm(distance_func, 4, lb, ub, options);

    % Output the optimal parameters and minimum distance
    u1_opt = uv_opt(1);
    v1_opt = uv_opt(2);
    u2_opt = uv_opt(3);
    v2_opt = uv_opt(4);
end

function distance = compute_surface_distance(params1, params2, u1, v1, u2, v2)
    % Get points on the surface of each superquadric
    [X1, Y1, Z1] = superquadric_point(params1, u1, v1);
    [X2, Y2, Z2] = superquadric_point(params2, u2, v2);
    
    % Calculate vector between the two points on the surfaces
    vec = [X1 - X2, Y1 - Y2, Z1 - Z2];
    
    % Calculate dot product to get squared distance
    distance = dot(vec, vec);
end





function [X, Y, Z] = superquadric_point(params, u, v)
    % Extract parameters for the superquadric
    a1 = params(1);
    a2 = params(2);
    a3 = params(3);
    epsilon1 = params(4);
    epsilon2 = params(5);
    theta = params(6);
    psi = params(7);
    phi = params(8);
    px = params(9);
    py = params(10);
    pz = params(11);

    % Transformation matrix T
    T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];

    % Parametric equations for the superquadric surface
    x = a1 * sign(cos(v)).*abs(cos(v)).^(epsilon2) .* sign(cos(u)).*abs(cos(u)).^(epsilon1);
    y = a2 * sign(cos(v)).*abs(cos(v)).^(epsilon2) .* sign(sin(u)).*abs(sin(u)).^(epsilon1);
    z = a3 * sign(sin(v)).*abs(sin(v)).^(epsilon2);

    % Ensure x, y, z are column vectors
    x = x(:);
    y = y(:);
    z = z(:);

    % Transform parametric coordinates
    superquadric_point = T * [x'; y'; z'; ones(size(x'))];

    % Extract X, Y, Z coordinates
    X = superquadric_point(1);
    Y = superquadric_point(2);
    Z = superquadric_point(3);
end



