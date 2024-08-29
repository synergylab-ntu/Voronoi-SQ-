function vertexGraph = Copy_of_createVoronoiGraph(polyhedrons, thresholdDistance, pointsPerFace)
    % Function to create a graph from a set of 3D polyhedrons, connect nodes based on proximity, and add uniformly distributed points on each face
    
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
        if size(vertices, 1) >= 4
            K = convhull(vertices);
        else
            K = [];
        end

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

        % Add uniformly distributed points on each face
        for k_idx = 1:size(K, 1)
            faceVertices = vertices(K(k_idx, :), :);
            
            % Generate grid points on the face
            facePoints = generateUniformPointsOnFace(faceVertices, pointsPerFace);
            
            % Add face points as nodes to the graph, ensuring uniqueness
            for fp_idx = 1:size(facePoints, 1)
                facePointStr = mat2str(facePoints(fp_idx, :));  % Convert face point to string to use as key
                if ~isKey(vertexIDMap, facePointStr)
                    vertexID = vertexGraph.numnodes + 1;
                    vertexIDMap(facePointStr) = vertexID;
                    vertexGraph = addnode(vertexGraph, table(vertexID, facePoints(fp_idx, 1), facePoints(fp_idx, 2), facePoints(fp_idx, 3), 'VariableNames', {'ID', 'X', 'Y', 'Z'}));
                end
            end
            
            % Add edges between the face points and the vertices on the face
            for fp_idx = 1:size(facePoints, 1)
                for v_idx = 1:size(faceVertices, 1)
                    startVertex = facePoints(fp_idx, :);
                    endVertex = faceVertices(v_idx, :);
                    
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
    end
    
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
    title('3D Graph with Uniform Face Points');
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

function facePoints = generateUniformPointsOnFace(faceVertices, pointsPerFace)
    % Function to generate uniformly distributed points on a triangular face
    
    % Barycentric coordinates for uniform sampling
    b = linspace(0, 1, pointsPerFace);
    [B1, B2] = meshgrid(b, b);
    valid = B1 + B2 <= 1;
    B1 = B1(valid);
    B2 = B2(valid);
    B3 = 1 - B1 - B2;
    
    % Calculate Cartesian coordinates of the points
    facePoints = B1 * faceVertices(1, :) + B2 * faceVertices(2, :) + B3 * faceVertices(3, :);
end

% function facePoints = generateUniformPointsOnFace(faceVertices, pointsPerFace)
%     % Function to generate uniformly distributed points on a triangular face
%     
%     % Vertices of the triangular face
%     v1 = faceVertices(1, :);
%     v2 = faceVertices(2, :);
%     v3 = faceVertices(3, :);
%     
%     % Generate random points within the triangle using random uniform sampling
%     u = rand(pointsPerFace, 1);
%     v = rand(pointsPerFace, 1);
%     
%     % Ensure points are within the triangle by adjusting coordinates
%     outside = u + v > 1;
%     u(outside) = 1 - u(outside);
%     v(outside) = 1 - v(outside);
%     
%     % Calculate Cartesian coordinates using barycentric coordinates
%     w = 1 - u - v;
%     facePoints = u .* v1 + v .* v2 + w .* v3;
% end










