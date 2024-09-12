function ShortestPath(vertexGraph, startPoint, endPoint, params)
    updatedGraph = vertexGraph;

    % Add temporary node at the closest point for the start point
    [startTempNodeID, startIntersection, updatedGraph] = addTemporaryNode3D(updatedGraph, startPoint);

    % Add temporary node at the closest point for the end point
    [endTempNodeID, endIntersection, updatedGraph] = addTemporaryNode3D(updatedGraph, endPoint);

    % Find the shortest path between the start and end temporary nodes
    pathNodes = shortestpath(updatedGraph, startTempNodeID, endTempNodeID);

    % Get coordinates of all nodes in the graph
    XData = updatedGraph.Nodes.X;
    YData = updatedGraph.Nodes.Y;
    ZData = updatedGraph.Nodes.Z;

    % Get coordinates of nodes along the shortest path
    XData_nodes = updatedGraph.Nodes.X(pathNodes);
    YData_nodes = updatedGraph.Nodes.Y(pathNodes);
    ZData_nodes = updatedGraph.Nodes.Z(pathNodes);

    % Prepend start point and append end point to path coordinates
    pathCoords = [startPoint; [XData_nodes, YData_nodes, ZData_nodes]; endPoint];

    % Plot the superquadrics
    figure;
    hold on;
    plot_multiple_superquadrics(params, 1);

    % Plot the shortest path
    plot3(pathCoords(:,1), pathCoords(:,2), pathCoords(:,3), 'r-', 'LineWidth', 2);

    % Plot nodes in the shortest path
    plot3(XData_nodes, YData_nodes, ZData_nodes, 'go', 'MarkerSize', 6,'MarkerFaceColor', 'b');

    % Plot start and end points, and intersection points
    scatter3(startPoint(1), startPoint(2), startPoint(3), 'g', 'filled');
    scatter3(endPoint(1), endPoint(2), endPoint(3), 'g', 'filled');
    scatter3(startIntersection(1), startIntersection(2), startIntersection(3), 'm', 'filled');
    scatter3(endIntersection(1), endIntersection(2), endIntersection(3), 'm', 'filled');

    % Plot lines from start/end points to the intersections on the edges
    plot3([startPoint(1), startIntersection(1)], [startPoint(2), startIntersection(2)], [startPoint(3), startIntersection(3)], 'k-', 'LineWidth', 2);
    plot3([endPoint(1), endIntersection(1)], [endPoint(2), endIntersection(2)], [endPoint(3), endIntersection(3)], 'k-', 'LineWidth', 2);

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    hold off;
end

function [tempNodeID, intersectionPoint, updatedGraph] = addTemporaryNode3D(vertexGraph, point)
    % Find the closest edge to the point and add a temporary node at the projection
    [startCoords, endCoords] = findClosestEdge3D(vertexGraph, point);
    
    % Compute vectors needed for projection
    edgeVector = endCoords - startCoords;
    pointVector = point - startCoords;
    
    % Project the point onto the edge
    edgeLengthSquared = sum(edgeVector.^2);
    t = dot(pointVector, edgeVector) / edgeLengthSquared;
    
    % Clamp t to the range [0, 1] to ensure the projection is within the edge segment
    t = max(0, min(1, t));
    
    % Compute the intersection point
    intersectionPoint = startCoords + t * edgeVector;
    
    % Generate a new node ID
    tempNodeID = vertexGraph.numnodes + 1;
    
    % Create a new table entry for the node
    newNodeTable = table(tempNodeID, intersectionPoint(1), intersectionPoint(2), intersectionPoint(3), ...
        'VariableNames', {'ID', 'X', 'Y', 'Z'});
    
    % Add the node to the graph
    updatedGraph = addnode(vertexGraph, newNodeTable);
    
    % Find the indices of the nodes closest to the intersection point
    startIdx = find(ismember([vertexGraph.Nodes.X, vertexGraph.Nodes.Y, vertexGraph.Nodes.Z], startCoords, 'rows'));
    endIdx = find(ismember([vertexGraph.Nodes.X, vertexGraph.Nodes.Y, vertexGraph.Nodes.Z], endCoords, 'rows'));
    
    % Compute distances for the new edges
    distStartToTemp = sqrt(sum((intersectionPoint - startCoords).^2));
    distTempToEnd = sqrt(sum((intersectionPoint - endCoords).^2));
    
    % Add new edges connecting the temporary node to the original edge's endpoints
    updatedGraph = addedge(updatedGraph, tempNodeID, startIdx, distStartToTemp);
    updatedGraph = addedge(updatedGraph, tempNodeID, endIdx, distTempToEnd);
end

function [startCoords, endCoords] = findClosestEdge3D(vertexGraph, point)
    % Initialize variables to track the closest edge
    minDistance = Inf;
    closestEdgeID = [];

    % Loop through each edge in the graph
    for edgeIdx = 1:numedges(vertexGraph)
        % Get the indices of the nodes at the ends of the edge
        startIdx = vertexGraph.Edges.EndNodes(edgeIdx, 1);
        endIdx = vertexGraph.Edges.EndNodes(edgeIdx, 2);
        
        % Get the coordinates of the start and end nodes
        startCoordsCurrent = [vertexGraph.Nodes.X(startIdx), vertexGraph.Nodes.Y(startIdx), vertexGraph.Nodes.Z(startIdx)];
        endCoordsCurrent = [vertexGraph.Nodes.X(endIdx), vertexGraph.Nodes.Y(endIdx), vertexGraph.Nodes.Z(endIdx)];
        
        % Calculate the squared distance from the point to the edge
        distSquared = pointToEdgeSquaredDistance3D(point, startCoordsCurrent, endCoordsCurrent);

        % Update the closest edge if the current edge is closer
        if distSquared < minDistance
            minDistance = distSquared;
            startCoords = startCoordsCurrent;
            endCoords = endCoordsCurrent;
        end
    end
end

function distSquared = pointToEdgeSquaredDistance3D(point, startCoords, endCoords)
    % Calculate the squared distance from a point to an edge in 3D
    x = point(1); y = point(2); z = point(3);
    x1 = startCoords(1); y1 = startCoords(2); z1 = startCoords(3);
    x2 = endCoords(1); y2 = endCoords(2); z2 = endCoords(3);

    % Compute vectors
    a = x - x1;
    b = y - y1;
    c = z - z1;
    d = x2 - x1;
    e = y2 - y1;
    f = z2 - z1;

    % Compute the squared length of the edge
    lenSq = d^2 + e^2 + f^2;

    % Calculate the projection parameter t
    if lenSq ~= 0 % Ensure the edge is not a point
        dotProd = a * d + b * e + c * f;
        t = dotProd / lenSq;
    else
        t = -1;
    end

    % Clamp t to the range [0, 1]
    if t < 0
        xx = x1; yy = y1; zz = z1;
    elseif t > 1
        xx = x2; yy = y2; zz = z2;
    else
        xx = x1 + t * d;
        yy = y1 + t * e;
        zz = z1 + t * f;
    end

    % Calculate the squared distance from the point to the closest point on the edge
    distSquared = (x - xx)^2 + (y - yy)^2 + (z - zz)^2;
end

