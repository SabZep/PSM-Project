
% create an empty structure with the desired fields
s = struct('name', [], 'xY', [], 'RT', [], 'dt', [], 'Y', [], 'psy', [], 'P', [], 'xn', [], 'ppi', []);

% loop through each PPI variable and concatenate data for each field
for i = 1:6
    % load the PPI struct from file
    filename = sprintf('PPI_PPI%d.mat', i);
    load(filename, 'PPI');

    % concatenate PPI data for each field
    s.name = [s.name; PPI.name];
    s.xY = [s.xY; PPI.xY];
    s.RT = [s.RT; PPI.RT];
    s.dt = [s.dt; PPI.dt];
    s.Y = [s.Y; PPI.Y];
    s.psy = [s.psy; PPI.psy];
    s.P = [s.P; PPI.P];
    s.xn = [s.xn; PPI.xn];
    s.ppi = [s.ppi; PPI.ppi];
end

% create a 1*1 struct with the concatenated data
result = struct('data', s);

% saving the data
save('concatenated_PPI.mat', 'result');