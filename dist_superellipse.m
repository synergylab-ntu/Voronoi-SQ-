t1 = linspace(-pi,pi,100);t2 = t1;d=zeros(100,100);
for i = 1:100
    for j = 1:100
        d(i,j) = distanceBetweenSuperellipses(t1(i), t2(j), params1, params2);
    end
end
surf(t1,t2,d')
xlabel('Angle parameter for first superellipse')
ylabel('Angle parameter for second superellipse')
zlabel('Euclidean Distance between points')

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