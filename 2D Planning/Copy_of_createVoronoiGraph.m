function vertexGraph = Copy_of_createVoronoiGraph(polyhedrons)
    % Function to create a graph from a set of polyhedrons and connect nodes based on proximity
    
    % Function to calculate Euclidean distance between two points
    euclideanDistance = @(p1, p2) sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2);

    % Step 1: Collect all vertices from polygons into a single matrix
    allVertices = [];
    
    % Loop through each polyhedron and concatenate their vertices into allVertices
    polygonsVertices = cell(1, size(polyhedrons, 1));
    for i = 1:size(polyhedrons, 1)
        rep = vrep(polyhedrons{i});
        vertices = unique((rep.V)', "rows");
        polygonsVertices{i} = orderVerticesCyclic(vertices);
        % Add the vertices to the overall list (each vertex as a row)
        allVertices = [allVertices; vertices];
    end
    
    % Step 2: Round the vertices to the 3rd decimal place
    allVertices = round(allVertices, 3);
    
    % Step 3: Remove duplicate vertices
    [uniqueVertices, ~] = unique(allVertices, 'rows');
    
    % Step 4: Initialize the vertex graph and ID map
    vertexIDMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    vertexGraph = graph();
    
    % Step 5: Add unique vertices as nodes to the graph
    numUniqueVertices = size(uniqueVertices, 1);
    for v_idx = 1:numUniqueVertices
        vertexStr = mat2str(uniqueVertices(v_idx, :));  % Convert vertex to string to use as key
        vertexID = vertexGraph.numnodes + 1;
        vertexIDMap(vertexStr) = vertexID;
        vertexGraph = addnode(vertexGraph, table(vertexID, uniqueVertices(v_idx, 1), uniqueVertices(v_idx, 2), 'VariableNames', {'ID', 'X', 'Y'}));
    end
    
    % Step 6: Initialize the map for tracking added edges
    addedEdges = containers.Map('KeyType', 'char', 'ValueType', 'logical');
    
    % Step 7: Add edges for each polygon
    for poly_idx = 1:length(polygonsVertices)
        vertices = polygonsVertices{poly_idx};
        numVertices = size(vertices, 1);
        
        % Add edges based on the polygon structure (connect consecutive vertices and close the loop)
        for e_idx = 1:numVertices
            % Round the start and end vertices to the 3rd decimal place
            startVertex = round(vertices(e_idx, :), 3);
            endVertex = round(vertices(mod(e_idx, numVertices) + 1, :), 3);  % Next vertex, wrapping around to close the loop
    
            % Convert them to strings for lookup in the map
            startStr = mat2str(startVertex);
            endStr = mat2str(endVertex);
    
            % Get the corresponding node IDs from vertexIDMap
            startIdx = vertexIDMap(startStr);
            endIdx = vertexIDMap(endStr);
    
            % Check for self-loops (i.e., startIdx ≠ endIdx)
            if startIdx ~= endIdx
                % Normalize the edge representation
                edgeStr = mat2str(sort([startIdx, endIdx]));
    
                % Calculate the Euclidean distance between the vertices
                distance = euclideanDistance(startVertex, endVertex);
    
                % Create a string that includes both the edge nodes and distance
                edgeKey = [edgeStr, '-', num2str(distance)];
    
                % Check if this edge with the same distance has already been added
                if ~isKey(addedEdges, edgeKey)
                    % Add the edge with the distance as the weight
                    vertexGraph = addedge(vertexGraph, startIdx, endIdx, distance);
    
                    % Mark this edge as added
                    addedEdges(edgeKey) = true;
                end
            end
        end
    end
    
    % Get all nodes from the graph
    allNodes = vertexGraph.Nodes;
    
    % Iterate over all pairs of nodes to check distances and add edges

    numNodes = height(allNodes);
    distances = [];

    % Step 2: Compute all pairwise distances between nodes
    for i = 1:numNodes
        for j = i+1:numNodes
            % Extract coordinates of node i and node j
            node_i_coords = [allNodes.X(i), allNodes.Y(i)];
            node_j_coords = [allNodes.X(j), allNodes.Y(j)];
            
            % Calculate Euclidean distance between node i and node j
            dist_ij = euclideanDistance(node_i_coords, node_j_coords);
            
            % Store the distances along with corresponding node indices
            distances = [distances; dist_ij, i, j]; % [distance, node i, node j]
        end
    end
    
    % Step 3: Sort the distances by the first column (the distances)
    distances = sortrows(distances, 1);

    thresholdDistance = min(vertexGraph.Edges.Weight);

    % Step 4: Add edges based on sorted distances and check for full connectivity
    for k = 1:size(distances, 1)
        node_i = distances(k, 2);
        node_j = distances(k, 3);
        dist_ij = distances(k, 1);

        % Check if the graph is fully connected
        if dist_ij >= thresholdDistance &&  max(conncomp(vertexGraph)) == 1 
            disp('Graph is fully connected.');
            break;
        end
        
        % Add the edge between the closest nodes
        vertexGraph = addedge(vertexGraph, node_i, node_j, dist_ij); 
    end

    % Plot the updated graph with edge labels (optional)
    figure;
    h = plot(vertexGraph, 'XData', vertexGraph.Nodes.X, 'YData', vertexGraph.Nodes.Y);
    h.NodeLabel = {};  % Remove node labels
    h.EdgeLabel = {};  % Remove edge labels
    title('Graph');
    xlabel('X');
    ylabel('Y');
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
