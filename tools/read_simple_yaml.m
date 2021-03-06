%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read_simple_yaml.m
%
% Extract from a simple and 1 level yaml the data in it
%
% Anthony Remazeilles and Jawad Masood
% Copyright Tecnalia 2019
% Beerware license.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = read_simple_yaml(filename)

    fid = fopen(filename);

    spec = "%s %s";
    infile = textscan(fid, spec, 'Delimiter', ':');

    labels = infile{1};
    values = infile{2};
    [nitem misc] = size(labels);

    for i=1:nitem
        label = strtrim(labels{i});
        if (length(label) ~= 0)
            data.(label) = values(i);
        end
    end
    fclose(fid);
end