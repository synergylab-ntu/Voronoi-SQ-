function Params = matrix_cell(superellipses)
% Convert to Params format
    Params = cell(size(superellipses, 1), 1);
    for i = 1:size(superellipses, 1)
        Params{i}.a = superellipses(i, 1);
        Params{i}.b = superellipses(i, 2);
        Params{i}.n = superellipses(i, 3);
        Params{i}.h = superellipses(i, 4);
        Params{i}.k = superellipses(i, 5);
        Params{i}.theta = superellipses(i, 6);
    end
end