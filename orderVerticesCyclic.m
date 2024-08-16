function orderedVertices = orderVerticesCyclic(vertices)
    % Compute convex hull
    K = convhull(vertices(:,1), vertices(:,2));
    
    % Extract ordered vertices from convex hull
    orderedVertices = vertices(K,:);
end