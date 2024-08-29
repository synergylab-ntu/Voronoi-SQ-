function vertexGraph = createVoronoiGraph(polyhedrons)
    % Function to create a graph from a set of 3D polyhedrons and connect nodes based on proximity
    
    % Function to calculate Euclidean distance between two points in 3D
    euclideanDistance = @(p1, p2) sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);
    
    % Initialize containers.Map to track unique vertices and their IDs
    vertexIDMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    vertexGraph = graph();
    
    % Initialize a set to track edges that have already been added
    edgeSet = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    
    % Add unique vertices as nodes to the graph
    for i = 1:size(polyhedrons, 1)
        rep = vrep(polyhedrons{i});
        vertices = unique((rep.V)', 'rows');
        
        % Calculate convex hull
        K = convhull(vertices);
        
        numVertices = size(vertices, 1);
        
        % Add vertices as nodes to the graph, ensuring uniqueness
        for v_idx = 1:numVertices
            vertexStr = mat2str(vertices(v_idx, :));  % Convert vertex to string to use as key
            if ~isKey(vertexIDMap, vertexStr)
                vertexID = vertexGraph.numnodes + 1;
                vertexIDMap(vertexStr) = vertexID;
                vertexGraph = addnode(vertexGraph, table(vertexID, vertices(v_idx, 1), vertices(v_idx, 2), vertices(v_idx, 3), 'VariableNames', {'ID', 'X', 'Y', 'Z'}));
            end
        end
        
        % Add edges based on the convex hull
        for k_idx = 1:size(K, 1)
            for col_idx = 1:size(K, 2)
                startVertex = vertices(K(k_idx, col_idx), :);
                endVertex = vertices(K(k_idx, mod(col_idx, size(K, 2)) + 1), :);
                
                startStr = mat2str(startVertex);
                endStr = mat2str(endVertex);
                
                startIdx = vertexIDMap(startStr);
                endIdx = vertexIDMap(endStr);
                
                % Create a sorted edge string to avoid duplicates
                edgeStr = mat2str(sort([startIdx, endIdx]));
                
                % Calculate the distance between the vertices
                distance = euclideanDistance(startVertex, endVertex);
                
                % Ensure start and end vertices are distinct before adding an edge
                if startIdx ~= endIdx && ~isKey(edgeSet, edgeStr)
                    vertexGraph = addedge(vertexGraph, startIdx, endIdx, distance);
                    edgeSet(edgeStr) = true;  % Mark the edge as added
                end
            end
        end
    end
    sortedWeights = sort(vertexGraph.Edges.Weight);
    thresholdDistance = 3*mean(sortedWeights(1:5));

    
    % Iterate over all pairs of nodes to check distances and add edges
    allNodes = vertexGraph.Nodes;
    numNodes = height(allNodes);
    for i = 1:numNodes
        for j = i+1:numNodes
            % Extract coordinates of node i and node j
            node_i_coords = [allNodes.X(i), allNodes.Y(i), allNodes.Z(i)];
            node_j_coords = [allNodes.X(j), allNodes.Y(j), allNodes.Z(j)];
            
            % Calculate Euclidean distance between node i and node j
            dist_ij = euclideanDistance(node_i_coords, node_j_coords);
            
            % Create a sorted edge string to avoid duplicates
            edgeStr = mat2str(sort([i, j]));
            
            % If distance is less than threshold and edge is not already added, add an edge between them
            if dist_ij <= thresholdDistance && ~isKey(edgeSet, edgeStr)
                vertexGraph = addedge(vertexGraph, i, j, dist_ij);
                edgeSet(edgeStr) = true;  % Mark the edge as added
            end
        end
    end
    
    % Plot the updated 3D graph (optional)
    figure;
    h = plot(vertexGraph, 'XData', vertexGraph.Nodes.X, 'YData', vertexGraph.Nodes.Y, 'ZData', vertexGraph.Nodes.Z, 'EdgeAlpha', 0.6);
    h.NodeLabel = {};  % Remove node labels
    h.EdgeLabel = {};  % Remove edge labels
    title('3D Graph');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    axis equal;

    % Check the connectivity of the graph using conncomp
    components = conncomp(vertexGraph);
    numComponents = max(components);
    if numComponents == 1
        disp('The graph is connected.');
    else
        disp(['The graph is not connected. It has ', num2str(numComponents), ' components.']);
    end
end




