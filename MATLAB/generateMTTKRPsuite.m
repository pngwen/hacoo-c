

% Generate a random sparse tensor, factor matrices, and mttkrp
% solutions over all modes.

% Set the filename to save the factor matrices
factorMatrixFilename = 'factor_matrices.txt';
mttkrpFilename = 'mttkrp_answers.txt';
sptensorFilename = 'sptensor_data.txt';

% Set number of columns
fmax = 3;

%generate random modes for sparse tensor
num_modes = 3;
MAX_MODE = 10; %max value for a mode

% Generate random modes
random_modes = zeros(1,num_modes);

for i = 1:num_modes
    random_modes(i) = randi([1, MAX_MODE]);
end

% Generate random sptensor
percent_nnz = 0.1; %percent of nnz in the tensor
tensor = generate_and_save_sptensor(random_modes, percent_nnz, sptensorFilename);

% Generate and save the factor matrices
U = generate_and_save_factor_matrices(random_modes, fmax, factorMatrixFilename);

% Calculate MTTKRP over all modes and save answers to file
mttkrp_ans = generate_and_save_mttkrp(tensor, U, mttkrpFilename);

% Function to generate random factor matrices and save them to a file
function U = generate_and_save_factor_matrices(tensor_dims, fmax, filename)
    % Input:
    % tensor_dims - A vector containing the dimensions of the tensor [I, J, K,...] 
    % fmax - factor matrices's max number of columns
    % filename - The name of the file to save the matrices

    % Example: tensor_dims = [I, J, K]; where I, J, K are the dimensions of the tensor

    % Number of modes (dimensions) of the tensor
    num_modes = length(tensor_dims);

    % Open file for writing
    fileID = open_file(filename);

    % Generate random factor matrices with dimensions compatible to the tensor
    U = cell(num_modes, 1);
    for i = 1:num_modes
        % Random matrix of size [tensor_dims(i), fmax] w/ vals b/w 0 and 10
        U{i} = 10 * rand(tensor_dims(i), fmax);

        % Write matrix dimensions before each matrix
        fprintf(fileID, '%d %d\n', tensor_dims(i), fmax);  % Write dimensions as [rows, columns]

        % Write the matrix values to the file, space-delimited
        [rows, cols] = size(U{i});
        for r = 1:rows
            fprintf(fileID, '%.6f ', U{i}(r, :));  % Write each row of the matrix
            fprintf(fileID, '\n');  % Move to the next line after each row
        end
    end

    % Close the file
    fclose(fileID);

    % Display success message
    fprintf('Random factor matrices have been generated and saved to "%s".\n', filename);
end

function mttkrp_ans = generate_and_save_mttkrp(tensor, U, filename)
    % Open file for writing
    fileID = open_file(filename);

    %to hold calculation answers
    mttkrp_ans = cell(ndims(tensor), 1);

    % Calculate MTTKRP over each mode & write to file
    for i = 1:ndims(tensor)
        mttkrp_ans{i} = mttkrp(tensor,U,i);
        %fprintf(fileID, 'MTTKRP for Mode %d:\n',i);

        %write dimensions for each matrix to file
        fprintf(fileID, '%d ', size( mttkrp_ans{i}));
        fprintf(fileID, '\n');

        [rows, cols] = size(mttkrp_ans{i});
        for r = 1:rows
            fprintf(fileID, '%.6f ', mttkrp_ans{i}(r, :));
            fprintf(fileID, '\n');
        end
    end

    % Close the file
    fclose(fileID);

    fprintf("Finished writing MTTKRP data to %s.\n", filename);
end

function tensor = generate_and_save_sptensor(modes, percent_nnz,filename)
    %Create a tensor with percent_nnz are nonzeros.
    tensor = sptenrand(modes,percent_nnz);

    %open and save tensor data to file
    fileID = open_file(filename);

    %write  modes first
    sz = size(tensor);
    for i=1:ndims(tensor)
        if i == ndims(tensor)
            %don't print space after last mode
            fprintf(fileID, '%d', sz(i));
        else
            fprintf(fileID, '%d ', sz(i));
        end
    end
    fprintf(fileID, '\n');

    idx_sz = size(tensor.subs);

    %write indexes and values
    for i = 1:idx_sz(1)
        fprintf(fileID, '%d ', tensor.subs(i, :)-1);
        fprintf(fileID, '%.6f', tensor.vals(i));
        fprintf(fileID, '\n');
    end

    % Close the file
    fclose(fileID);

    fprintf("Finished writing tensor data to %s.\n", filename);
end

%open file and return fileID
function fileID = open_file(filename) 
    % Open file for writing
    fileID = fopen(filename, 'w');  % Open the file in write mode

    if fileID == -1
        error('Cannot open the file %s for writing.', filename);
    end
end