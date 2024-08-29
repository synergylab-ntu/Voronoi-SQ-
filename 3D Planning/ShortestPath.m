function ShortestPath(vertexGraph, startPoint, endPoint,params)
    % Find closest nodes to the start and end points
    startNode = findClosestNode(vertexGraph, startPoint);
    endNode = findClosestNode(vertexGraph, endPoint);
    
    % Find shortest path between the closest nodes
    pathNodes = shortestpath(vertexGraph, startNode, endNode);
    
    % Plot the graph with highlighted shortest path and start/end lines
    figure;
    h = plot(vertexGraph, 'XData', vertexGraph.Nodes.X, 'YData', vertexGraph.Nodes.Y, 'ZData', vertexGraph.Nodes.Z);
    h.NodeLabel = {};  % Remove node labels
    h.EdgeLabel = {};  % Remove edge labels
    title('3D Graph with Shortest Path Highlighted');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;
    hold on;
    
    % Highlight the shortest path
    highlight(h, pathNodes, 'EdgeColor', 'r', 'LineWidth', 2);
    
    % Plot lines from start and end points to their closest nodes
    plot3([startPoint(1), vertexGraph.Nodes.X(startNode)], [startPoint(2), vertexGraph.Nodes.Y(startNode)], [startPoint(3), vertexGraph.Nodes.Z(startNode)], 'k--', 'LineWidth', 1.5);
    plot3([endPoint(1), vertexGraph.Nodes.X(endNode)], [endPoint(2), vertexGraph.Nodes.Y(endNode)], [endPoint(3), vertexGraph.Nodes.Z(endNode)], 'k--', 'LineWidth', 1.5);
    
    % Mark start and end points
    plot3(startPoint(1), startPoint(2), startPoint(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    plot3(endPoint(1), endPoint(2), endPoint(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    
    %plot superquadrics
    plot_multiple_superquadrics(params);
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

