clear all;close all;
superellipses = [
    0.5, 0.75, 1.2, 1.5, 1, pi/3;
    0.7, 0.5, 1, -1, 1, pi/6;
    2, 1, 0.2, -2.5, -2, -pi/6;
    0.6, 0.8, 0.7, 1, -1, pi/5;
    0.5, 0.5, 1, -0.6, -0.65, 0];
params = matrix_cell(superellipses);
[x_grid, y_grid] = meshgrid(linspace(-5, 5, 200), linspace(-5, 5, 200));
%V_overall = se_contour_overall(x_grid, y_grid, params);
%surf(x_grid,y_grid,V_overall')
n_se = size(superellipses,1);
params_robot = matrix_cell([0.3,0.7,0.2, -1, 0,pi/2]);

% Find optimal orientation 
params_robot.theta = find_optimal_orientation(params_robot, params);

% Generate superellipse contours for visualization
t = linspace(-pi, pi, 100);
for i = 1:n_se
    paramsi = params{i};
    [xp, yp] = superellipse(t, paramsi.a, paramsi.b, paramsi.n);
    [Xp(i,:), Yp(i,:)] = transform(xp, yp, paramsi.theta, paramsi.h, paramsi.k);
end

[xr, yr] = superellipse(t, params_robot.a, params_robot.b, params_robot.n);
[Xr(:), Yr(:)] = transform(xr, yr, params_robot.theta, params_robot.h, params_robot.k);

axis equal
hold on

for i = 1:n_se
    fill(Xp(i,:), Yp(i,:), 'b', 'FaceAlpha', 0.3); % Adjust FaceAlpha for transparency
end
fill(Xr(:), Yr(:), 'r', 'FaceAlpha', 0.3);
hold off

function V_overall = se_contour_overall(x, y, params)
    % Initialize the overall potential as 0 (additive identity)
    V_overall = 0;

    % Loop through each superellipse and compute the sum of potentials
    for i = 1:length(params)
        F = in_out(x, y, params{i});
        V = 10^6 * exp(-0.5 * F);
        V_overall = V_overall + V;
    end
end

function F = in_out(x,y,params)
    a = params.a;b=params.b;n=params.n;h=params.h;k=params.k;theta = params.theta;
    % Reverse translation
    x_trans = x - h;
    y_trans = y - k;

    % Reverse rotation by angle -theta
    x_rot = x_trans * cos(-theta) - y_trans * sin(-theta);
    y_rot = x_trans * sin(-theta) + y_trans * cos(-theta);

    % Superellipse implicit function
    F = (abs(x_rot/a).^(2/n)) + (abs(y_rot/b).^(2/n)) - 1;
end


function [Fx_overall, Fy_overall] = tot_grad(x, y, params)
    % Initialize the overall gradient components
    Fx_overall = 0;
    Fy_overall = 0;
    
    % Loop through each superellipse to compute the gradient contributions
    for i = 1:length(params)
        % Get individual parameters
        a = params{i}.a;
        b = params{i}.b;
        n = params{i}.n;
        h = params{i}.h;
        k = params{i}.k;
        theta = params{i}.theta;
        
        % Reverse translation
        x_trans = x - h;
        y_trans = y - k;
        
        % Reverse rotation by angle -theta
        x_rot = x_trans * cos(-theta) - y_trans * sin(-theta);
        y_rot = x_trans * sin(-theta) + y_trans * cos(-theta);
        
        % Superellipse implicit function
        F = (abs(x_rot/a).^(2/n)) + (abs(y_rot/b).^(2/n)) - 1;
        
        % Compute the gradient of F with respect to x_rot and y_rot
        dF_dx_rot = (2/n) * (abs(x_rot/a).^(2/n - 1)) .* sign(x_rot) / a;
        dF_dy_rot = (2/n) * (abs(y_rot/b).^(2/n - 1)) .* sign(y_rot) / b;
        
        % Chain rule to get the gradient with respect to x and y
        dF_dx = dF_dx_rot * cos(-theta) + dF_dy_rot * sin(-theta);
        dF_dy = -dF_dx_rot * sin(-theta) + dF_dy_rot * cos(-theta);
        
        % Compute the gradient of the individual potential function
        V = 10^4 * exp(-2 * F);
        Vx = -2 * 10^4 * exp(-2 * F) .* dF_dx;
        Vy = -2 * 10^4 * exp(-2 * F) .* dF_dy;
        
        % Sum the gradients of individual potentials to get the overall gradient
        Fx_overall = Fx_overall + Vx;
        Fy_overall = Fy_overall + Vy;
    end
end

function resultant_moment = compute_resultant_moment(params_robot, Params)
    t = linspace(-pi,pi,100);
    % Get points on the boundary of the robot
    [xr, yr] = superellipse(t, params_robot.a, params_robot.b, params_robot.n);
    [Xr, Yr] = transform(xr, yr, params_robot.theta, params_robot.h, params_robot.k);

    M_total = 0;
    center_x = params_robot.h;  % Reference point (center of the robot)
    center_y = params_robot.k;

    % Sum the moments acting on each boundary point
    for i = 1:length(Xr)
        [Fx, Fy] = tot_grad(Xr(i), Yr(i), Params);
        % Position vector from the robot's center to the boundary point
        rx = Xr(i) - center_x;
        ry = Yr(i) - center_y;

        % Compute the moment (torque) as the cross product of r and F
        M = rx * Fy - ry * Fx;
        M_total = M_total + M;
    end

    % Resultant moment
    resultant_moment = M_total;
end


function optimal_orientation = find_optimal_orientation(params_robot, Params)
    orientations = linspace(0, 2*pi, 100); % Check multiple orientations
    min_moment = inf;
    optimal_orientation = params_robot.theta;

    for i = 1:length(orientations)
        params_robot.theta = orientations(i);
        moment = compute_resultant_moment(params_robot, Params);

        if abs(moment) < min_moment
            min_moment = abs(moment);
            optimal_orientation = orientations(i);
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

% function M_theta = resultant_moment_theta(theta, params_robot, Params)
%     % Update the robot's orientation
%     params_robot.theta = theta;
%     % Compute the resultant moment for this orientation
%     M_theta = compute_resultant_moment(params_robot, Params);
% end
% 
% function optimal_orientation = find_optimal_orientation(params_robot, Params)
%     % Define the function handle for resultant moment as a function of theta
%     moment_fun = @(theta) abs(resultant_moment_theta(theta, params_robot, Params));
% 
%     % Define the optimization options
%     options = optimoptions('particleswarm', 'Display', 'iter', ...
%                            'SwarmSize', 50, 'MaxIterations', 100);
% 
%     % Use particleswarm to find the optimal orientation (minimizing the absolute moment)
%     [optimal_orientation] = particleswarm(moment_fun, 1, 0, 2*pi, options);
% end



