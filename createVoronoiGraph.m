function vertexGraph = createVoronoiGraph(polyhedrons, thresholdDistance)
    % Function to create a graph from a set of polyhedrons and connect nodes based on proximity
    
    % Function to calculate Euclidean distance between two points
    euclideanDistance = @(p1, p2) sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2);
    
    % Generate the vertices for each polyhedron and store them in polygonsVertices
    polygonsVertices = cell(1, size(polyhedrons, 1));
    for i = 1:size(polyhedrons, 1)
        rep = vrep(polyhedrons{i});
        vertices = unique((rep.V)', "rows");
        polygonsVertices{i} = orderVerticesCyclic(vertices);
    end
    
    % Initialize containers.Map to track unique vertices and their IDs
    vertexIDMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    vertexGraph = graph();
    
    % Add unique vertices as nodes to the graph
    for poly_idx = 1:length(polygonsVertices)
        vertices = polygonsVertices{poly_idx};
        numVertices = size(vertices, 1);
        
        % Add vertices as nodes to the graph, ensuring uniqueness
        for v_idx = 1:numVertices
            vertexStr = mat2str(vertices(v_idx, :));  % Convert vertex to string to use as key
            if ~isKey(vertexIDMap, vertexStr)
                vertexID = vertexGraph.numnodes + 1;
                vertexIDMap(vertexStr) = vertexID;
                vertexGraph = addnode(vertexGraph, table(vertexID, vertices(v_idx, 1), vertices(v_idx, 2), 'VariableNames', {'ID', 'X', 'Y'}));
            end
        end
        
        % Add edges based on the polygon structure (connect consecutive vertices and close the loop)
        for e_idx = 1:numVertices
            startVertex = vertices(e_idx, :);
            endVertex = vertices(mod(e_idx, numVertices) + 1, :);  % Next vertex, wrapping around for closing the loop
            
            startStr = mat2str(startVertex);
            endStr = mat2str(endVertex);
            
            startIdx = vertexIDMap(startStr);
            endIdx = vertexIDMap(endStr);
            
            % Calculate the distance between the vertices
            distance = euclideanDistance(startVertex, endVertex);
            
            % Ensure start and end vertices are distinct before adding an edge
            if startIdx ~= endIdx
                vertexGraph = addedge(vertexGraph, startIdx, endIdx, distance);
            end
        end
    end
    
    
    % Get all nodes from the graph
    allNodes = vertexGraph.Nodes;
    
    % Iterate over all pairs of nodes to check distances and add edges
    numNodes = height(allNodes);
    for i = 1:numNodes
        for j = i+1:numNodes
            % Extract coordinates of node i and node j
            node_i_coords = [allNodes.X(i), allNodes.Y(i)];
            node_j_coords = [allNodes.X(j), allNodes.Y(j)];
            
            % Calculate Euclidean distance between node i and node j
            dist_ij = euclideanDistance(node_i_coords, node_j_coords);
            
            % If distance is less than threshold, add an edge between them
            if dist_ij <= thresholdDistance
                vertexGraph = addedge(vertexGraph, i, j, dist_ij);
            end
        end
    end

    
    % Display message about adding edges
    %disp(['Edges added between nodes within ', num2str(thresholdDistance), ' units of each other.']);
    
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
