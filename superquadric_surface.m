function superquadric_surface(a1, a2, a3, epsilon1, epsilon2, theta, psi, phi, px, py, pz)
    % Create the transformation matrix T
    T = [
        cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
        sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
        -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
        0, 0, 0, 1
    ];
    
    % Extract n, o, a vectors from T
    n = T(1:3, 1);
    o = T(1:3, 2);
    a = T(1:3, 3);
    p = T(1:3, 4);
    
    % Define the parametric equations
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

    % Plot the superquadric surface
    fsurf(@(u,v) parametric_eq(u,v), [-pi/2, pi/2, -pi, pi], 'MeshDensity', 50);
    
    axis equal;
    title('Superquadric Surface');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    camlight;
    lighting gouraud;
    
    % Adjust the view to see more of the shape
    view(30, 30);
    
    % Add a grid for better perspective
    grid on;
end