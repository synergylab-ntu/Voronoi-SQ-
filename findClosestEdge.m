function closestEdgeID = findClosestEdge(vertexGraph, point)
    minDistance = Inf;
    closestEdgeID = [];
    for edgeIdx = 1:numedges(vertexGraph)
        % Get the indices of the nodes at the ends of the edge
        startIdx = vertexGraph.Edges.EndNodes(edgeIdx, 1);
        endIdx = vertexGraph.Edges.EndNodes(edgeIdx, 2);
        
        % Get the coordinates of the start and end nodes
        startCoords = [vertexGraph.Nodes.X(startIdx), vertexGraph.Nodes.Y(startIdx)];
        endCoords = [vertexGraph.Nodes.X(endIdx), vertexGraph.Nodes.Y(endIdx)];
        
        % Calculate the squared distance from the point to the edge
        distSquared = pointToEdgeSquaredDistance(point, startCoords, endCoords);
        
        % Update the closest edge if the current edge is closer
        if distSquared < minDistance
            minDistance = distSquared;
            closestEdgeID = edgeIdx;
        end
    end
end

function distSquared = pointToEdgeSquaredDistance(point, startCoords, endCoords)
    edgeVector = endCoords - startCoords;
    pointVector = point - startCoords;
    edgeLengthSquared = sum(edgeVector.^2);
    
    if edgeLengthSquared == 0
        % The edge is actually a point
        distSquared = sum((point - startCoords).^2);
        return;
    end
    
    % Calculate the projection parameter t
    t = dot(pointVector, edgeVector) / edgeLengthSquared;
    t = max(0, min(1, t));
    
    % Calculate the projection of the point onto the edge
    projection = startCoords + t * edgeVector;
    
    % Calculate the squared distance from the point to the projection
    distSquared = sum((point - projection).^2);
end
