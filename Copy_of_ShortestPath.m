function Copy_of_ShortestPath(vertexGraph, startPoint, endPoint, params)
    % Find closest nodes to the start and end points
    startNode = findClosestNode(vertexGraph, startPoint);
    endNode = findClosestNode(vertexGraph, endPoint);
    
    % Find shortest path between the closest nodes
    pathNodes = shortestpath(vertexGraph, startNode, endNode);
    
    % Extract coordinates of nodes in the shortest path
    pathX = vertexGraph.Nodes.X(pathNodes);
    pathY = vertexGraph.Nodes.Y(pathNodes);
    pathZ = vertexGraph.Nodes.Z(pathNodes);
    
    % Plot the superquadrics
    figure;
    hold on;
    plot_multiple_superquadrics(params,1);
    
    % Plot lines from start and end points to their closest nodes
    plot3([startPoint(1), vertexGraph.Nodes.X(startNode)], [startPoint(2), vertexGraph.Nodes.Y(startNode)], [startPoint(3), vertexGraph.Nodes.Z(startNode)], 'r-', 'LineWidth', 2);
    plot3([endPoint(1), vertexGraph.Nodes.X(endNode)], [endPoint(2), vertexGraph.Nodes.Y(endNode)], [endPoint(3), vertexGraph.Nodes.Z(endNode)], 'r-', 'LineWidth', 2);
    
    % Mark start and end points
    plot3(startPoint(1), startPoint(2), startPoint(3), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    plot3(endPoint(1), endPoint(2), endPoint(3), 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
    
    % Highlight the shortest path
    plot3(pathX, pathY, pathZ, 'r-', 'LineWidth', 2);
    
    % Plot nodes in the shortest path
    plot3(pathX, pathY, pathZ, 'go', 'MarkerSize', 6,'MarkerFaceColor', 'b');
    
    % Set axis labels and equal scaling
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    
    % Adjust view angle if needed
    view(3); % Set to 3D view
end

function closestNode = findClosestNode(vertexGraph, point)
    % Function to find the closest node in vertexGraph to a given 3D point
    
    % Extract node coordinates
    nodeCoords = [vertexGraph.Nodes.X, vertexGraph.Nodes.Y, vertexGraph.Nodes.Z];
    
    % Calculate distances to all nodes using pdist2
    distances = pdist2(point, nodeCoords);
    
    % Find the index of the node with the minimum distance
    [~, minIdx] = min(distances);
    
    % Return the closest node ID
    closestNode = minIdx;
end

