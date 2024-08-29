function plotClosestPoints3D(params1, params2)
    % Find the closest points between two superquadrics
    [xc1, yc1, zc1, xc2, yc2, zc2] = findClosestPoints3D(params1, params2);
    
    % Plot the superquadrics and closest points
    plotSuperquadrics(params1);
    hold on;
    plotSuperquadrics(params2);
    plot3(xc1, yc1, zc1, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    plot3(xc2, yc2, zc2, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
    plot3([xc1, xc2], [yc1, yc2], [zc1, zc2], 'k--', 'LineWidth', 2);
    hold off;
    
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Closest Points between Two Superquadrics');
    axis equal;
    grid on;
end

function plotSuperquadrics(params)
    % Extract parameters
    a1 = params{1}(1);
    a2 = params{1}(2);
    a3 = params{1}(3);
    e1 = params{1}(4);
    e2 = params{1}(5);
    theta = params{1}(6);
    psi = params{1}(7);
    phi = params{1}(8);
    px = params{1}(9);
    py = params{1}(10);
    pz = params{1}(11);
    
    % Plot the superquadric surface
    superquadric_surface(a1, a2, a3, theta, psi, phi, px, py, pz, e1, e2);
end

function [xc1, yc1, zc1, xc2, yc2, zc2] = findClosestPoints3D(params1, params2)
    % Objective function for optimization
    function d = distanceBetweenSuperquadrics(t)
        % Extract angles
        t1 = t(1:2);
        t2 = t(3:4);
        
        % Calculate points on the superquadrics
        [x1, y1, z1] = superquadric(t1, params1{1}(1), params1{1}(2), params1{1}(3), params1{1}(4), params1{1}(5));
        [x2, y2, z2] = superquadric(t2, params2{1}(1), params2{1}(2), params2{1}(3), params2{1}(4), params2{1}(5));
        
        % Transform points according to superquadric parameters
        [xc1, yc1, zc1] = transform3D(x1, y1, z1, params1{1}(6), params1{1}(7), params1{1}(8), params1{1}(9), params1{1}(10), params1{1}(11));
        [xc2, yc2, zc2] = transform3D(x2, y2, z2, params2{1}(6), params2{1}(7), params2{1}(8), params2{1}(9), params2{1}(10), params2{1}(11));
        
        % Calculate the Euclidean distance
        d = sqrt((xc2 - xc1)^2 + (yc2 - yc1)^2 + (zc2 - zc1)^2);
    end

    % Extract centers of superquadrics for initial guess calculation
    h1 = params1{1}(9); h2 = params2{1}(9);
    k1 = params1{1}(10); k2 = params2{1}(10);
    l1 = params1{1}(11); l2 = params2{1}(11);
    dz = l2 - l1; dy = k2 - k1; dx = h2 - h1;
    
    % Initial guess
    x0 = [atan2(dy, dx), atan2(dz, sqrt(dx^2 + dy^2)), atan2(-dy, -dx), atan2(-dz, -sqrt(dx^2 + dy^2))];
    
    % Lower and upper bounds for u and v
    lb = [-pi/2, -pi,-pi/2,-pi];  % Lower bounds for u and v
    ub = [pi/2, pi,pi/2,pi];    % Upper bounds for u and v

    % Optimization options for particleswarm
    options = optimoptions('particleswarm', 'Display', 'off', 'UseParallel', false);%, 'InitialSwarmMatrix', x0, 'SwarmSize', 100);

    % Run the global optimization using particleswarm
    [t_opt, ~] = particleswarm(@distanceBetweenSuperquadrics, 4, lb, ub, options);
    %[t_opt, ~] = fmincon(@distanceBetweenSuperquadrics, x0, [], [], [], [], lb, ub, [], options);

%     % Calculate the closest points using the optimized parameters
%     [x1, y1, z1] = superquadric(t_opt(1:2), params1.a, params1.b, params1.c, params1.e1, params1.e2);
%     [x2, y2, z2] = superquadric(t_opt(3:4), params2.a, params2.b, params2.c, params2.e1, params2.e2);
%     [xc1, yc1, zc1] = transform3D(x1, y1, z1, params1.theta, params1.phi, params1.psi, params1.h, params1.k, params1.l);
%     [xc2, yc2, zc2] = transform3D(x2, y2, z2, params2.theta, params2.phi, params2.psi, params2.h, params2.k, params2.l);
     % Calculate points on the superquadrics
        [x1, y1, z1] = superquadric(t_opt(1:2), params1{1}(1), params1{1}(2), params1{1}(3), params1{1}(4), params1{1}(5));
        [x2, y2, z2] = superquadric(t_opt(3:4), params2{1}(1), params2{1}(2), params2{1}(3), params2{1}(4), params2{1}(5));
        
        % Transform points according to superquadric parameters
        [xc1, yc1, zc1] = transform3D(x1, y1, z1, params1{1}(6), params1{1}(7), params1{1}(8), params1{1}(9), params1{1}(10), params1{1}(11));
        [xc2, yc2, zc2] = transform3D(x2, y2, z2, params2{1}(6), params2{1}(7), params2{1}(8), params2{1}(9), params2{1}(10), params2{1}(11));
end

function [x, y, z] = superquadric(t, a, b, c, e1, e2)
    % Superquadric surface parameterization
    u = t(1);
    v = t(2);
    x = a * sign(cos(v)) * abs(cos(v))^e2 * sign(cos(u)) * abs(cos(u))^e1;
    y = b * sign(cos(v)) * abs(cos(v))^e2 * sign(sin(u)) * abs(sin(u))^e1;
    z = c * sign(sin(v)) * abs(sin(v))^e2;
end

function [X, Y, Z] = transform3D(x, y, z, theta, phi, psi, px, py, pz)
T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];
%     % Apply rotation
%     Rz = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
%     Ry = [cos(phi), 0, sin(phi); 0, 1, 0; -sin(phi), 0, cos(phi)];
%     Rx = [1, 0, 0; 0, cos(psi), -sin(psi); 0, sin(psi), cos(psi)];
%     R = Rz * Ry * Rx;
%     p = R * [x; y; z];
%     
%     % Apply translation
%     x = p(1) + h;
%     y = p(2) + k;
%     z = p(3) + l;
 % Ensure x, y, z are column vectors
    x = x(:);
    y = y(:);
    z = z(:);

    % Transform parametric coordinates
    superquadric_point = T * [x'; y'; z'; ones(size(x'))];

    % Extract X, Y, Z coordinates
    X = reshape(superquadric_point(1,:), size(x));
    Y = reshape(superquadric_point(2,:), size(x));
    Z = reshape(superquadric_point(3,:), size(x));

end

function superquadric_surface(a1, a2, a3, theta, psi, phi, px, py, pz, e1, e2)
    % Create the transformation matrix T
    T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];
    
    % Define the superquadric function
    function p_w = F(u, v)
        x_w = a1 * sign(cos(u)) * abs(cos(u))^e2 * sign(cos(v)) * abs(cos(v))^e1;
        y_w = a2 * sign(cos(u)) * abs(cos(u))^e2 * sign(sin(v)) * abs(sin(v))^e1;
        z_w = a3 * sign(sin(u)) * abs(sin(u))^(2/2);
        
        % Apply transformation
        p_w = T * [x_w; y_w; z_w; 1];
    end
    
    % Create grid
    [U, V] = meshgrid(linspace(-pi, pi, 100), linspace(-pi/2, pi/2, 100));
    
    % Calculate surface
    X = zeros(size(U));
    Y = zeros(size(U));
    Z = zeros(size(U));
    for i = 1:numel(U)
        P = F(U(i), V(i));
        X(i) = P(1);
        Y(i) = P(2);
        Z(i) = P(3);
    end
    
    % Plot surface
    surf(X, Y, Z, 'FaceAlpha', 0.5);
    shading interp;
end

