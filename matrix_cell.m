function Params = matrix_cell(superellipses)
    if size(superellipses, 1) == 1
        % If superellipses is a single row, return a struct
        Params.a = superellipses(1, 1);
        Params.b = superellipses(1, 2);
        Params.n = superellipses(1, 3);
        Params.h = superellipses(1, 4);
        Params.k = superellipses(1, 5);
        Params.theta = superellipses(1, 6);
    else
        % If superellipses has multiple rows, return a cell array of structs
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
end
