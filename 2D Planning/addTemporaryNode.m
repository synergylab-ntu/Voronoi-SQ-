function [tempNodeID, intersectionPoint, updatedGraph] = addTemporaryNode(vertexGraph, point)
    % Use findClosestEdge to get the closest edge's start and end coordinates
    [startCoords, endCoords] = findClosestEdge(vertexGraph, point);
    
    % Compute the vectors needed for projection
    edgeVector = endCoords - startCoords;
    pointVector = point - startCoords;
    
    % Calculate the parameter t for the projection of the point onto the edge
    edgeLengthSquared = sum(edgeVector.^2);
    t = dot(pointVector, edgeVector) / edgeLengthSquared;
    
    % Clamp t to the range [0, 1]
    t = max(0, min(1, t));
    
    % Calculate the intersection point
    intersectionPoint = startCoords + t * edgeVector;
    
    % Generate a new node ID
    tempNodeID = vertexGraph.numnodes + 1;
    
    % Create a new table entry for the new node
    newNodeTable = table(tempNodeID, intersectionPoint(1), intersectionPoint(2), 'VariableNames', {'ID', 'X', 'Y'});
    
    % Add the new node to the graph
    updatedGraph = addnode(vertexGraph, newNodeTable);
    
    % Get the indices of the nodes closest to the intersection point
    startIdx = find(ismember([vertexGraph.Nodes.X, vertexGraph.Nodes.Y], startCoords, 'rows'));
    endIdx = find(ismember([vertexGraph.Nodes.X, vertexGraph.Nodes.Y], endCoords, 'rows'));
    
    % Calculate the distances (weights) for the new edges
    distStartToTemp = sqrt(sum((intersectionPoint - startCoords).^2));
    distTempToEnd = sqrt(sum((intersectionPoint - endCoords).^2));
    
    % Add new edges connecting the new node to the original edge's endpoints
    updatedGraph = addedge(updatedGraph, tempNodeID, startIdx, distStartToTemp);
    updatedGraph = addedge(updatedGraph, tempNodeID, endIdx, distTempToEnd);
    
    % Optionally, remove the original edge if you no longer need it
    % updatedGraph = rmedge(updatedGraph, edgeIdx);

end
