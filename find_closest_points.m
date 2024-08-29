% function [u1_opt, v1_opt, u2_opt, v2_opt, min_distance] = find_closest_points(params1, params2)
%     % Define the distance function to be minimized
%     distance_func = @(uv) norm(superquadric_point(params1, uv(1), uv(2)) - superquadric_point(params2, uv(3), uv(4)));
% 
%     % Initial guess for u and v within bounds
%     initial_guess = [0, 0, 0, 0];
% 
%     % Lower and upper bounds for u and v
%     lb = [-pi/2, -pi, -pi/2, -pi];
%     ub = [pi/2, pi, pi/2, pi];
% 
%     % Optimization options for particleswarm
%     options = optimoptions('particleswarm', 'Display', 'off', 'UseParallel', false, 'SwarmSize', 100);
% 
%     % Perform optimization to find optimal u and v
%     [uv_opt, min_distance] = particleswarm(distance_func, 4, lb, ub, options);
% 
%     % Output the optimal parameters and minimum distance
%     u1_opt = uv_opt(1);
%     v1_opt = uv_opt(2);
%     u2_opt = uv_opt(3);
%     v2_opt = uv_opt(4);
% 
%     % Plotting the superquadrics and closest points
%     plot_superquadrics(params1, params2, u1_opt, v1_opt, u2_opt, v2_opt);
% end
% 
% function plot_superquadrics(params1, params2, u1_opt, v1_opt, u2_opt, v2_opt)
%     % Generate grid of u and v values for plotting
%     [U, V] = meshgrid(linspace(-pi/2, pi/2, 50), linspace(-pi, pi, 50));
% 
%     % Plot superquadric 1
%     [X1, Y1, Z1] = generate_superquadric_surface(params1, U, V);
%     surf(X1, Y1, Z1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%     hold on;
% 
%     % Plot superquadric 2
%     [X2, Y2, Z2] = generate_superquadric_surface(params2, U, V);
%     surf(X2, Y2, Z2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
% 
%     % Calculate and plot closest points
%     P1_opt = superquadric_point(params1, u1_opt, v1_opt);
%     P2_opt = superquadric_point(params2, u2_opt, v2_opt);
%     plot3([P1_opt(1), P2_opt(1)], [P1_opt(2), P2_opt(2)], [P1_opt(3), P2_opt(3)], 'ro-', 'LineWidth', 2, 'MarkerSize', 10);
% 
%     % Additional plot settings
%     axis equal;
%     title('Superquadric Surfaces and Closest Points');
%     xlabel('X'); ylabel('Y'); zlabel('Z');
%     % legend('Superquadric 1', 'Superquadric 2', 'Closest Points');
%     view(3);
%     grid on;
%     camlight;
%     lighting gouraud;
%     hold off;
% end
% 
% function [X, Y, Z] = generate_superquadric_surface(params, U, V)
%     % Allocate space for the coordinates
%     X = zeros(size(U));
%     Y = zeros(size(U));
%     Z = zeros(size(U));
% 
%     % Compute the coordinates for each point on the grid
%     for i = 1:numel(U)
%         P = superquadric_point(params, U(i), V(i));
%         X(i) = P(1);
%         Y(i) = P(2);
%         Z(i) = P(3);
%     end
% end
% 
% function P = superquadric_point(params, u, v)
%     % Extract parameters for the superquadric
%     a1 = params{1}(1);
%     a2 = params{1}(2);
%     a3 = params{1}(3);
%     epsilon1 = params{1}(4);
%     epsilon2 = params{1}(5);
%     theta = params{1}(6);
%     psi = params{1}(7);
%     phi = params{1}(8);
%     px = params{1}(9);
%     py = params{1}(10);
%     pz = params{1}(11);
% 
%     % Parametric equations for the superquadric surface
%     x = a1 * sign(cos(v)) .* abs(cos(v)).^epsilon2 .* sign(cos(u)) .* abs(cos(u)).^epsilon1;
%     y = a2 * sign(cos(v)) .* abs(cos(v)).^epsilon2 .* sign(sin(u)) .* abs(sin(u)).^epsilon1;
%     z = a3 * sign(sin(v)) .* abs(sin(v)).^epsilon2;
% 
%     % Transformation matrix T
%     T = [
%         cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
%         sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
%         -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
%         0, 0, 0, 1
%     ];
% 
%     % Apply the transformation
%     P_homogeneous = T * [x; y; z; 1];
%     P = P_homogeneous(1:3);
% end

function [u1_opt, v1_opt, u2_opt, v2_opt, min_distance, plane_eq] = find_closest_points(params1, params2)
    % Define the distance function to be minimized
    distance_func = @(uv) norm(superquadric_point(params1, uv(1), uv(2)) - superquadric_point(params2, uv(3), uv(4)));

    % Initial guess for u and v within bounds
    initial_guess = [0, 0, 0, 0];

    % Lower and upper bounds for u and v
    lb = [-pi/2, -pi, -pi/2, -pi];
    ub = [pi/2, pi, pi/2, pi];

    % Optimization options for particleswarm
    options = optimoptions('particleswarm', 'Display', 'off', 'UseParallel', false, 'SwarmSize', 100);

    % Perform optimization to find optimal u and v
    [uv_opt, min_distance] = particleswarm(distance_func, 4, lb, ub, options);

    % Output the optimal parameters and minimum distance
    u1_opt = uv_opt(1);
    v1_opt = uv_opt(2);
    u2_opt = uv_opt(3);
    v2_opt = uv_opt(4);

    % Calculate the closest points
    P1_opt = superquadric_point(params1, u1_opt, v1_opt);
    P2_opt = superquadric_point(params2, u2_opt, v2_opt);
    
    % Calculate the midpoint of the closest points
    midpoint = (P1_opt + P2_opt) / 2;
    
    % Calculate the normal vector of the plane
    normal_vector = P2_opt - P1_opt;

    % Calculate the plane equation coefficients
    d = -dot(normal_vector, midpoint);
    plane_eq = [normal_vector; d];

    % Plotting the superquadrics, closest points, and bisecting plane
    plot_superquadrics(params1, params2, u1_opt, v1_opt, u2_opt, v2_opt, plane_eq);
end

function plot_superquadrics(params1, params2, u1_opt, v1_opt, u2_opt, v2_opt, plane_eq)
    % Generate grid of u and v values for plotting
    [U, V] = meshgrid(linspace(-pi/2, pi/2, 50), linspace(-pi, pi, 50));

    % Plot superquadric 1
    [X1, Y1, Z1] = generate_superquadric_surface(params1, U, V);
    surf(X1, Y1, Z1, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    hold on;

    % Plot superquadric 2
    [X2, Y2, Z2] = generate_superquadric_surface(params2, U, V);
    surf(X2, Y2, Z2, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    % Calculate and plot closest points
    P1_opt = superquadric_point(params1, u1_opt, v1_opt);
    P2_opt = superquadric_point(params2, u2_opt, v2_opt);
    plot3([P1_opt(1), P2_opt(1)], [P1_opt(2), P2_opt(2)], [P1_opt(3), P2_opt(3)], 'ro-', 'LineWidth', 2, 'MarkerSize', 10);

    % Plot the bisecting plane
    plot_bisecting_plane(plane_eq);

    % Additional plot settings
    axis equal;
    axis([-2 6 -2.5 4 -2 7])
    title('Superquadric Surfaces, Closest Points, and Bisecting Plane');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    % legend('Superquadric 1', 'Superquadric 2', 'Closest Points', 'Bisecting Plane');
    view(3);
    grid on;
    camlight;
    lighting gouraud;
    hold off;
end

function plot_bisecting_plane(plane_eq)
    [x, y] = meshgrid(linspace(-5, 7, 10), linspace(-5, 5, 10));
    z = (-plane_eq(1) * x - plane_eq(2) * y - plane_eq(4)) / plane_eq(3);
    surf(x, y, z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', 'g');
end

function [X, Y, Z] = generate_superquadric_surface(params, U, V)
    % Allocate space for the coordinates
    X = zeros(size(U));
    Y = zeros(size(U));
    Z = zeros(size(U));

    % Compute the coordinates for each point on the grid
    for i = 1:numel(U)
        P = superquadric_point(params, U(i), V(i));
        X(i) = P(1);
        Y(i) = P(2);
        Z(i) = P(3);
    end
end

function P = superquadric_point(params, u, v)
    % Extract parameters for the superquadric
    a1 = params{1}(1);
    a2 = params{1}(2);
    a3 = params{1}(3);
    epsilon1 = params{1}(4);
    epsilon2 = params{1}(5);
    theta = params{1}(6);
    psi = params{1}(7);
    phi = params{1}(8);
    px = params{1}(9);
    py = params{1}(10);
    pz = params{1}(11);

    % Parametric equations for the superquadric surface
    x = a1 * sign(cos(v)) .* abs(cos(v)).^epsilon2 .* sign(cos(u)) .* abs(cos(u)).^epsilon1;
    y = a2 * sign(cos(v)) .* abs(cos(v)).^epsilon2 .* sign(sin(u)) .* abs(sin(u)).^epsilon1;
    z = a3 * sign(sin(v)) .* abs(sin(v)).^epsilon2;

    % Transformation matrix T
    T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];

    % Apply the transformation
    P_homogeneous = T * [x; y; z; 1];
    P = P_homogeneous(1:3);
end





