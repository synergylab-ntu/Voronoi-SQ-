function boundaryPoints = segmentImageAndExtractBoundaries(imagePath,sam)
    % Load the image
    image = imread(imagePath);

    % Display the image and let the user select foreground points
    imshow(image);
    title('Click on the foreground points, then press Enter');
    [foregroundX, foregroundY] = ginput;
    
    % Store the foreground points in a matrix
    foregroundPoints = [foregroundX, foregroundY];
    close(gcf); % Close the image window after input

    % Segment the image using SAM
    imageSize = size(image);
    embeddings = extractEmbeddings(sam, image); % Extract feature embeddings from the image
    masks = segmentObjectsFromEmbeddings(sam, embeddings, imageSize, ...
        ForegroundPoints=foregroundPoints); % Create binary mask specifying object occupancy

    % Extract boundary points for each object in the mask
    boundaryPoints = bwboundaries(masks);

    % Determine the number of points in each boundary
    numPoints = cellfun(@(x) size(x, 1), boundaryPoints);

    % Sort the boundaries based on the number of points in descending order
    [~, sortIdx] = sort(numPoints, 'descend');

    % Keep only the top m boundaries
    m = size(foregroundPoints, 1); % Define the number of boundaries you want to keep
    if m > length(sortIdx)
        m = length(sortIdx); % Adjust m if it's larger than the number of available boundaries
    end

    % Select the top m boundaries
    boundaryPoints = boundaryPoints(sortIdx(1:m));
    %Transforming boundary points from image frame to cartesian frame
    for k = 1:length(boundaryPoints)
        boundary = boundaryPoints{k}; % Get the k-th boundary
        boundaryPoints{k} = [boundary(:,2), -boundary(:,1)];
    end
end