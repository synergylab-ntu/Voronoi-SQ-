% Function to add start and end nodes perpendicularly from start and end points to closest edges,
% find shortest path between them, and visualize
function pathCoords = ShortestPath(vertexGraph, startPoint, endPoint,Params)
    % Find closest edge indices
    startEdgeIdx = findClosestEdge(vertexGraph, startPoint);
    endEdgeIdx = findClosestEdge(vertexGraph, endPoint);
    
    % Add temporary node at perpendicular reach for start point
    [startTempNodeID, startIntersection, updatedGraph] = addTemporaryNode(vertexGraph, startEdgeIdx, startPoint);
    
    % Add temporary node at perpendicular reach for end point
    [endTempNodeID, endIntersection, updatedGraph] = addTemporaryNode(updatedGraph, endEdgeIdx, endPoint);
    
    % Find shortest path between start and end temporary nodes
    pathNodes = shortestpath(updatedGraph, startTempNodeID, endTempNodeID);
    
    % Get coordinates of nodes for plotting
    XData = updatedGraph.Nodes.X;
    YData = updatedGraph.Nodes.Y;

    % Get coordinates of nodes for the shortest path
    XData_nodes = updatedGraph.Nodes.X(pathNodes);
    YData_nodes = updatedGraph.Nodes.Y(pathNodes);
    
    % Prepend start point and append end point
    pathCoords = [startPoint; [XData_nodes, YData_nodes]; endPoint];

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
    
    % Plot the graph with specified node positions and edges on the existing figure
    hold on;
    h = plot(updatedGraph, 'XData', XData, 'YData', YData, 'EdgeAlpha', 0.6);
     % Hide node labels
    h.NodeLabel = {};
    
    % Plot start and end points, temporary nodes, and perpendicular paths
    scatter(startPoint(1), startPoint(2), 'g', 'filled');
    scatter(endPoint(1), endPoint(2), 'g', 'filled');
    scatter(startIntersection(1), startIntersection(2), 'm', 'filled');
    scatter(endIntersection(1), endIntersection(2), 'm', 'filled');
    plot([startPoint(1), startIntersection(1)], [startPoint(2), startIntersection(2)], 'k-', 'LineWidth', 2);
    plot([endPoint(1), endIntersection(1)], [endPoint(2), endIntersection(2)], 'k-', 'LineWidth', 2);

    for i = 1:n_se
         fill(Xp(i,:), Yp(i,:), 'b', 'FaceAlpha', 0.3); % Adjust FaceAlpha for transparency
    end
    
    % Highlight shortest path between start and end nodes
    highlight(h, pathNodes, 'EdgeColor', 'k', 'LineWidth', 2);
    
    % Customize plot appearance
    title('Shortest Path');
   % legend('Graph', 'Start Point', 'End Point', 'Start Node', 'End Node', 'Location', 'Best');
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
