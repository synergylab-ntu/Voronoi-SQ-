function vertexGraph = createVoronoiGraph(polyhedrons, initialThresholdDistance, distanceStep)
    % Function to create a graph from a set of 3D polyhedrons and connect nodes based on proximity

    % Function to calculate Euclidean distance between two points in 3D
    euclideanDistance = @(p1, p2) sqrt((p1(1) - p2(1))^2 + (p1(2) - p2(2))^2 + (p1(3) - p2(3))^2);

    % Initialize containers.Map to track unique vertices and their IDs
    vertexIDMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    vertexGraph = graph();

    % Initialize a set to track edges that have already been added
    edgeSet = containers.Map('KeyType', 'char', 'ValueType', 'logical');

    % Add unique vertices as nodes to the graph
    allVertices = [];
    for i = 1:size(polyhedrons, 1)
        rep = vrep(polyhedrons{i});
        vertices = unique((rep.V)', 'rows');
        allVertices = [allVertices; vertices]; % Collect all vertices

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

    % Function to compute minimum threshold distance needed to connect the graph fully using binary search
    function minThresholdDistance = computeMinThresholdDistanceBinarySearch(vertices)
        % Determine the minimum and maximum possible distances
        numVertices = size(vertices, 1);
        minDistance = inf;
        maxDistance = 0;
        
        for i = 1:numVertices
            for j = i+1:numVertices
                dist = euclideanDistance(vertices(i, :), vertices(j, :));
                if dist < minDistance
                    minDistance = dist;
                end
                if dist > maxDistance
                    maxDistance = dist;
                end
            end
        end
        
        % Binary search for the minimum threshold distance
        while maxDistance - minDistance > 1e-3  % Tolerance level
            midDistance = (minDistance + maxDistance) / 2;
            tempGraph = vertexGraph;

            % Add edges based on midDistance
            for i = 1:numVertices
                for j = i+1:numVertices
                    dist = euclideanDistance(vertices(i, :), vertices(j, :));
                    if dist <= midDistance
                        startIdx = vertexIDMap(mat2str(vertices(i, :)));
                        endIdx = vertexIDMap(mat2str(vertices(j, :)));
                        edgeStr = mat2str(sort([startIdx, endIdx]));
                        if ~isKey(edgeSet, edgeStr)
                            tempGraph = addedge(tempGraph, startIdx, endIdx, dist);
                            edgeSet(edgeStr) = true;
                        end
                    end
                end
            end

            % Check the connectivity of the temporary graph using conncomp
            components = conncomp(tempGraph);
            numComponents = max(components);
            if numComponents == 1
                maxDistance = midDistance;  % Decrease upper bound
            else
                minDistance = midDistance;  % Increase lower bound
            end
        end
        
        minThresholdDistance = maxDistance;  % The smallest distance that connects the graph
    end

    % Compute the minimum threshold distance needed to connect the graph fully
    minThresholdDistance = computeMinThresholdDistanceBinarySearch(allVertices);

    % Iterate over all pairs of nodes to check distances and add edges based on the minimum threshold distance
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

            % If distance is less than or equal to minThresholdDistance and edge is not already added, add an edge between them
            if dist_ij <= minThresholdDistance && ~isKey(edgeSet, edgeStr)
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
    
    % Display the minimum threshold distance needed to connect the graph
    disp(['The minimum threshold distance needed to connect the graph is ', num2str(minThresholdDistance)]);
end
