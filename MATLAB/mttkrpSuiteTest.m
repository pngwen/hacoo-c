%open file and return fileID
function fileID = open_file(filename) 
    % Open file for writing
    fileID = fopen(filename, 'w');  % Open the file in write mode

    if fileID == -1
        error('Cannot open the file %s for writing.', filename);
    end
end

% Set the filename to save the factor matrices
factorMatrixFilename = 'factor_matrices.txt';
mttkrpFilename = 'mttkrp_answers.txt';
sptensorFilename = 'sptensor_data.txt';

%------------------------Set up tensor
subs = [0,1,0;
        1,0,1;
        0,2,1]; %<-- Subscripts of the nonzeros.
subs = subs + 1;
vals = [1; 2; 3]; %<-- The values of the nonzeros.

tensor = sptensor(subs,vals); %<--the tensor

%write sptensor data to file
fileID = open_file(sptensorFilename);

%write modes to file
idx_sz = size(tensor.subs);
fprintf(fileID, '%d ', size(tensor));
fprintf(fileID, '\n');

%write indexes and values
for i = 1:idx_sz(1)
    fprintf(fileID, '%d ', tensor.subs(i, :));
    fprintf(fileID, '%.2f', tensor.vals(i));
    fprintf(fileID, '\n');
end


% Close the file
fclose(fileID);

%------------------------Set up factor matrices
fmax = 3;

a = [1,3,5;2,4,6];
b = [1,4,7;2,5,8;3,6,9];
c = [1,2,3;4,5,6];

U = {a,b,c}; %<--the cell array

%Write factor matrices to file
fileID = open_file(factorMatrixFilename);

    for i = 1:ndims(tensor)
        % Write matrix dimensions before each matrix
        fprintf(fileID, '%d %d\n', size(U{i},1), fmax);  % Write dimensions as [rows, columns]

        % Write the matrix values to the file, space-delimited
        [rows, cols] = size(U{i});
        for r = 1:rows
            fprintf(fileID, '%d ', U{i}(r, :));  % Write each row of the matrix
            fprintf(fileID, '\n');  % Move to the next line after each row
        end
    end

% Close the file
fclose(fileID);

%------------------------Write mttrkrp results to file
% Open file
fileID = open_file(mttkrpFilename);

%to hold calculation answers
mttkrp_ans = cell(ndims(tensor), 1);

% Calculate MTTKRP over each mode & write to file
for i = 1:ndims(tensor)
    mttkrp_ans{i} = mttkrp(tensor,U,i);
    %fprintf(fileID, 'MTTKRP for Mode %d:\n',i);
    [rows, cols] = size(mttkrp_ans{i});

    %write dimensions for each matrix
    fprintf(fileID, '%d ', size( mttkrp_ans{i}));
    fprintf(fileID, '\n');

    for r = 1:rows
        %write matrix rows
        fprintf(fileID, '%.2f ', mttkrp_ans{i}(r, :));
        fprintf(fileID, '\n');
    end
end

% Close the file
fclose(fileID);
