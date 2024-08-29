function plot_multiple_superquadrics(params,method)
    % Main function body
    %figure;
    hold on;

    for i = 1:length(params)
        % Extract parameters for the current superquadric
        a1 = params{i}(1);
        a2 = params{i}(2);
        a3 = params{i}(3);
        epsilon1 = params{i}(4);
        epsilon2 = params{i}(5);
        theta = params{i}(6);
        psi = params{i}(7);
        phi = params{i}(8);
        px = params{i}(9);
        py = params{i}(10);
        pz = params{i}(11);

        % Define transformation matrix T
        T = [
            cos(phi)*cos(theta)*cos(psi) - sin(phi)*sin(psi), -cos(phi)*cos(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*sin(theta), px;
            sin(phi)*cos(theta)*cos(psi) + cos(phi)*sin(psi), -sin(phi)*cos(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*sin(theta), py;
            -sin(theta)*cos(psi), sin(theta)*sin(psi), cos(theta), pz;
            0, 0, 0, 1
        ];

        % Parametric representation of superquadric surface
        u = linspace(-pi/2, pi/2, 50);  % Adjust the number of points as needed
        v = linspace(-pi, pi, 50);      % Adjust the number of points as needed
        [U, V] = meshgrid(u, v);

        x = a1 * sign(cos(V)).*abs(cos(V)).^(epsilon2) .* sign(cos(U)).*abs(cos(U)).^(epsilon1);
        y = a2 * sign(cos(V)).*abs(cos(V)).^(epsilon2) .* sign(sin(U)).*abs(sin(U)).^(epsilon1);
        z = a3 * sign(sin(V)).*abs(sin(V)).^(epsilon2);

        % Transform parametric coordinates
        X = T(1,1)*x + T(1,2)*y + T(1,3)*z + T(1,4);
        Y = T(2,1)*x + T(2,2)*y + T(2,3)*z + T(2,4);
        Z = T(3,1)*x + T(3,2)*y + T(3,3)*z + T(3,4);

        % Plot the superquadric surface
        if method == 0 % fill
            surf(X, Y, Z, 'EdgeColor', 'none');
        elseif method == 1 % do not fill
            surf(X,Y,Z, 'EdgeColor', 'none','FaceAlpha', 0.5)
        end
        % Optional: Plot mesh lines if needed
        % mesh(U, V, reshape([X;Y;Z], [], length(U)));
    end

    axis equal;
    title('Multiple Superquadric Surfaces');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    view(3);
    camlight;
    lighting gouraud;
    %hold off;
end
