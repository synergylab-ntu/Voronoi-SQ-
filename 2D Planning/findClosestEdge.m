function [startCoords, endCoords] = findClosestEdge(vertexGraph, point)
    % Initialize variables to track the closest edge
    minDistance = Inf;
    closestEdgeID = [];
    
    % Loop through each edge in the graph
    for edgeIdx = 1:numedges(vertexGraph)
        % Get the indices of the nodes at the ends of the edge
        startIdx = vertexGraph.Edges.EndNodes(edgeIdx, 1);
        endIdx = vertexGraph.Edges.EndNodes(edgeIdx, 2);
        
        % Get the coordinates of the start and end nodes
        startCoordsCurrent = [vertexGraph.Nodes.X(startIdx), vertexGraph.Nodes.Y(startIdx)];
        endCoordsCurrent = [vertexGraph.Nodes.X(endIdx), vertexGraph.Nodes.Y(endIdx)];
        
        % Calculate the squared distance from the point to the edge
        distSquared = pointToEdgeSquaredDistance(point, startCoordsCurrent, endCoordsCurrent);

        % Update the closest edge if the current edge is closer
        if distSquared < minDistance
            minDistance = distSquared;
            closestEdgeID = edgeIdx;
            startCoords = startCoordsCurrent;
            endCoords = endCoordsCurrent;
        end
    end
end

function distSquared = pointToEdgeSquaredDistance(point, startCoords, endCoords)
    % Unpack point and coordinates
    x = point(1);
    y = point(2);
    x1 = startCoords(1);
    y1 = startCoords(2);
    x2 = endCoords(1);
    y2 = endCoords(2);
    
    % Calculate vectors
    a = x - x1;
    b = y - y1;
    c = x2 - x1;
    d = y2 - y1;
    
    % Calculate squared length of the segment
    lenSq = c^2 + d^2;
    
    % Calculate the parameter t for projection
    if lenSq ~= 0 % In case of zero length line segment
        dot = a * c + b * d;
        t = dot / lenSq;
    else
        t = -1;
    end
    
    % Find the closest point on the line segment
    if t < 0
        xx = x1;
        yy = y1;
    elseif t > 1
        xx = x2;
        yy = y2;
    else
        xx = x1 + t * c;
        yy = y1 + t * d;
    end
    
    % Calculate the distance from the point to the closest point on the line segment
    dx = x - xx;
    dy = y - yy;
    distSquared = dx^2 + dy^2;
end
