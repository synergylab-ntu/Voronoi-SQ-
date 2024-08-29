function vertexGraph =  connectGraphWithMST(vertexGraph)
    % Compute Minimum Spanning Tree (MST)
    mstEdges = minspantree(vertexGraph);

    % Add MST edges to the original vertexGraph
    for i = 1:numedges(mstEdges)
        edgeNodes = mstEdges.Edges.EndNodes(i,:);
        node1 = edgeNodes(1);
        node2 = edgeNodes(2);

        % Retrieve the weight (distance) of the MST edge
        weight = mstEdges.Edges.Weight(i);

        % Add edge to vertexGraph with the calculated weight
        vertexGraph = addedge(vertexGraph, node1, node2, weight);
    end
end



