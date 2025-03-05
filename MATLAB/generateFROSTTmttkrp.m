% Load sparse tensor in COO format, generate random factor matrices, 
% write mttkrp solutions over all modes.
%Converts indices to 0-based indices before calculating MTTKRP 

% Set the filename to save the factor matrices
factorMatrixFilename = 'factor_matrices.txt';
mttkrpFilename = 'mttkrp_answers.txt';

% Set number of columns
fmax = 4;

[subs, vals] = read_coo("chicago.tns");

tensor = sptensor(subs,vals);

% Generate and save the factor matrices
U = generate_and_save_factor_matrices(size(tensor), fmax, factorMatrixFilename);

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
        mttkrp_ans{i} = mttkrp_0_based(tensor,U,i);
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

function M = mttkrp_0_based(X, U, mode)
    % X: Sparse tensor in COO format
    % U: Cell array of factor matrices {U1, U2, ..., Un}
    % mode: Mode along which to compute MTTKRP
    
    % Extract COO representation
    subs = X.subs - 1; % Convert to 0-based indices
    vals = X.vals;
    R = size(U{1}, 2); % Rank of factor matrices
    N = ndims(X);      % Number of modes
    
    % Initialize MTTKRP result matrix
    M = zeros(size(U{mode}, 1), R);
    
    % Iterate over non-zero elements
    for i = 1:length(vals)
        idx = subs(i, :); % 0-based indices
        temp = vals(i);
        
        for n = 1:N
            if n ~= mode
                temp = temp .* U{n}(idx(n) + 1, :); % Convert back to 1-based
            end
        end
        
        % Accumulate into the output matrix
        M(idx(mode) + 1, :) = M(idx(mode) + 1, :) + temp;
    end
end

%open file and return fileID
function fileID = open_file(filename) 
    % Open file for writing
    fileID = fopen(filename, 'w');  % Open the file in write mode

    if fileID == -1
        error('Cannot open the file %s for writing.', filename);
    end
end

% Read a COO text file & return indexes & values as separate arrays
function [idx, vals] = read_coo(filename)
    file = filename;

    %Get the first line using fgetl to figure out how many modes
    opt = {'Delimiter',' '};
    fid = fopen(file,'rt');
    hdr = fgetl(fid);
    num = numel(regexp(hdr,' ','split'));
    if strcmp(file,"enron.txt") || strcmp(file,"nell.txt") || strcmp(file,"lbnl.txt") || strcmp(file,"nell-2.txt") || strcmp(file,"lbnl-network.txt")
        fmt = repmat('%d',1,num-1); %to read files with decimal values (enron, nell-2,lbnl)
        fmt = strcat(fmt,'%f');
    else
        fmt = repmat('%d',1,num); %to read files with no decimal values
    end

    frewind(fid); %put first line back

    sizeA = [num Inf];
    tdata = fscanf(fid,fmt,sizeA);
    tdata = tdata';

    fclose(fid);

    idx = tdata(:,1:num-1);
    vals = tdata(:,end);
end
