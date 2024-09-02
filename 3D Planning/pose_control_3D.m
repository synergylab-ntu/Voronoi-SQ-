close all;clear;
%Define parameters for 11 superquadrics (a,b,c,e1,e2,phi,psi,theta,h,k,l)
params = {
    [1.25, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0];
    [2, 1, 2, 1, 1, pi/4, pi/4, pi/4, 4, 2, 2];
    [1, 2, 1, 0.5, 0.5, pi/3, pi/3, pi/3, -9, 0, -2];
    [1, 1, 2, 0.3, 0.4, pi/6, pi/6, pi/6, 5, -5, 5];
    [2, 3, 2, 1.5, 1.5, 0, pi/2, 0, -6, 7, -6];
    [1, 2.5, 1.5, 0.3, 0.5, pi/4, 0, pi/4, 5, 6.5, 3];
    [1.5, 1, 2, 1.5, 1.5, pi/3, pi/3, pi/3, -4, 2, 0];
    [2, 1.5, 3, 1, 0.5, pi/6, pi/6, pi/6, -9, -6, 3];
    [1.5, 1, 1.5, 0.5, 1, 0, pi/4, 0, 1, -4, -4];
    [1, 2, 2.5, 1.5, 0.3, pi/4, 0, pi/4, 0, 5, -2];
    [1.5, 2.5, 2, 1, 2, pi/3, pi/3, pi/3, -3, -4, 7]
};
params_robot = [1,3,2,0.9,1,0,0,0,11,12,10];

function F = in_out(x, y, z, params)
    % params : (a,b,c,e1,e2,phi,psi,theta,h,k,l)
    % (x, y, z) Point where F is to be calculated
    
    % Extract parameters
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

    % Rotation matrix R based on Euler angles (phi, theta, psi)
    R = [cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta);
         sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta);
         -sin(theta)*cos(psi),                              sin(theta)*sin(psi),                               cos(theta)];

    % Translation vector
    p = [px; py; pz];

    % Inverse of the rotation matrix (transpose of R since R is orthogonal)
    R_inv = R';

    % Inverse of the translation vector
    p_inv = -R_inv * p;

    % Construct the inverse transformation matrix T_inv
    T_inv = [R_inv, p_inv;
             0, 0, 0, 1];
    val = T_inv * [x; y; z; 1];
    x_loc = val(1);
    y_loc = val(2);
    z_loc = val(3);

    % Define the modified power function to handle signs
    pow = @(base, exp) (abs(base).^exp) .* sign(base);

    % Compute F using the custom power function
    F = pow(pow(x_loc/a1, 2/epsilon2) + pow(y_loc/a2, 2/epsilon2), epsilon2/epsilon1) + pow(z_loc/a3, 2/epsilon1) - 1;
end

function neg_grad_V = neg_grad_V(x, y, z, params)
    % params : (a,b,c,e1,e2,phi,psi,theta,h,k,l)
    % (x, y, z) Point where the negative gradient of V is to be calculated

    % Extract parameters
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

    % Rotation matrix R based on Euler angles (phi, theta, psi)
    R = [cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta);
         sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta);
         -sin(theta)*cos(psi),                              sin(theta)*sin(psi),                               cos(theta)];

    % Translation vector
    p = [px; py; pz];

    % Inverse of the rotation matrix (transpose of R since R is orthogonal)
    R_inv = R';

    % Inverse of the translation vector
    p_inv = -R_inv * p;

    % Construct the inverse transformation matrix T_inv
    T_inv = [R_inv, p_inv;
            0, 0, 0, 1];
    val = T_inv*[x;y;z;1];
    x_loc = val(1); y_loc = val(2); z_loc = val(3);

    % Calculate the partial derivatives of F with respect to x_loc, y_loc, z_loc
    F_x = (2/epsilon1) .* Pow(x_loc./a1, (2/epsilon2 - 1))./ (a1 .* epsilon2) .* ...
          Pow(Pow(x_loc./a1, 2/epsilon2) + Pow(y_loc./a2, 2/epsilon2), (epsilon2/epsilon1 - 1));

    F_y = (2/epsilon1) .* Pow(y_loc./a2, (2/epsilon2 - 1)) / (a2 * epsilon2) .* ...
          Pow(Pow(x_loc./a1, 2/epsilon2) + Pow(y_loc./a2, 2/epsilon2), (epsilon2/epsilon1 - 1));

    F_z = (2/epsilon1) .* Pow(z_loc./a3, (2/epsilon1 - 1)) / a3;

    % Gradient of F in local coordinates
    grad_F_local = [F_x; F_y; F_z];

    % Transform gradient to global coordinates
    grad_F_global = R * grad_F_local;

    % Compute the potential function V = 10^4 * exp(-2F)
    F = Pow(Pow(x_loc/a1, 2/epsilon2) + Pow(y_loc/a2, 2/epsilon2), epsilon2/epsilon1) + ...
        Pow(z_loc/a3, 2/epsilon1) - 1;
    V = 10^4 * exp(-2 * F);

    % Compute the gradient of V using the chain rule: grad(V) = -2V * grad(F)
    grad_V = -2 .* V .* grad_F_global;

    % Return the negative gradient
    neg_grad_V = -grad_V;
end


function result = Pow(base, exp)
    if exp == floor(exp)  % Check if exp is an integer
        % If exp is an integer, perform standard power operation
        result = base.^exp;
    else
        % If exp is non-integer, perform signed power operation
        result = (abs(base).^exp) .* sign(base);
    end
end

function [Fx, Fy, Fz] = tot_grad(x, y, z, params)
    % Calculate the total force based on the sum of the potential fields from multiple superquadrics
    % Input:
    %   x, y, z   : Coordinates where the force is to be calculated
    %   params    : Cell array where each cell contains parameters for one superquadric
    % Output:
    %   Fx, Fy, Fz : Components of the total force at the given coordinates
    
    % Initialize the total gradient components to zero
    Fx_total = 0;
    Fy_total = 0;
    Fz_total = 0;
    
    % Iterate over each set of parameters in the cell array
    for i = 1:length(params)
        % Get the current set of parameters
        param_set = params{i};
        
        % Compute the negative gradient of the potential field for the current superquadric
        neg_grad = neg_grad_V(x, y, z, param_set);

        % Accumulate the gradients
        Fx_total = Fx_total + neg_grad(1);
        Fy_total = Fy_total + neg_grad(2);
        Fz_total = Fz_total + neg_grad(3);
    end
    
    % Output the total force components
    Fx = Fx_total;
    Fy = Fy_total;
    Fz = Fz_total;
end

function resultant_moment = compute_resultant_moment(params_robot, params_superquadrics)
    % params_robot: Parameters of the robot superquadric
    % params_superquadrics: Cell array of parameters of all other superquadrics

    % Extract robot parameters
    a1 = params_robot(1);
    a2 = params_robot(2);
    a3 = params_robot(3);
    epsilon1 = params_robot(4);
    epsilon2 = params_robot(5);
    theta = params_robot(6);
    psi = params_robot(7);
    phi = params_robot(8);
    px = params_robot(9);
    py = params_robot(10);
    pz = params_robot(11);
    
    % Transformation matrix for robot
    T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];

    % Define the parametric equations for the superquadric surface
    function [x, y, z] = parametric_eq(u, v)
        x_w = a1 * sign(cos(u)) .* abs(cos(u)).^(2/epsilon2) .* sign(cos(v)) .* abs(cos(v)).^(2/epsilon1);
        y_w = a2 * sign(cos(u)) .* abs(cos(u)).^(2/epsilon2) .* sign(sin(v)) .* abs(sin(v)).^(2/epsilon1);
        z_w = a3 * sign(sin(u)) .* abs(sin(u)).^(2/epsilon2);
        
        % Apply transformation
        transformed = T * [x_w; y_w; z_w; ones(size(x_w))];
        x = transformed(1, :);
        y = transformed(2, :);
        z = transformed(3, :);
    end

    % Define ranges for u and v
    u_range = linspace(-pi/2, pi/2, 100);  % Sampling 100 points for u
    v_range = linspace(-pi, pi, 100);      % Sampling 100 points for v

    % Initialize total moment
    resultant_moment = [0; 0; 0];
    
    % Iterate over the surface points
    for u = u_range
        for v = v_range
            % Get the surface point
            [x, y, z] = parametric_eq(u, v);
            point = [x; y; z];
            
            % Calculate the total force at this point due to all superquadrics
            [Fx, Fy, Fz] = tot_grad(x, y, z, params_superquadrics); 
            
            % Compute the position vector relative to the center of the robot
            r = point - [px; py; pz];
            
            % Compute the moment at this point: M = r x F
            moment = cross(r, [Fx;Fy;Fz]);
            
            % Sum the moment contributions
            resultant_moment = resultant_moment + moment;
        end
    end
end
