clc; clear;
setPath();

tic
filename = 'D:\Jahandar\Lab\images\immune\results\MantonBM1\outs\filtered_gene_bc_matrices_h5.h5';
cutoff_thres = 100;       % remove genes that the sum of counts is below thres
dim_red_method = 'SNE';   % Dimension Reduction Method 'PCA' | 'SNE'
red_dim = 3;              % reduce the dataset features to this dim
normalize = 0;

%% Load data
fprintf('=============================================================\n');
fprintf('Reading the data...\n');
fprintf('input file = %s\n', filename);
hinfo = hdf5info(filename);

data.barcodes = hdf5read(filename, '/GRCh38/barcodes');
data.data = hdf5read(filename, '/GRCh38/data');
data.gene_names = hdf5read(filename, '/GRCh38/gene_names');
data.genes = hdf5read(filename, '/GRCh38/genes');
data.indices = hdf5read(filename, '/GRCh38/indices') + 1;
data.indptr = hdf5read(filename, '/GRCh38/indptr') + 1;
data.shape = hdf5read(filename, '/GRCh38/shape');

% change gene_names to array
gene_names = cell(length(data.gene_names), 1);
for i = 1:length(data.gene_names)
    gene_names{i} = data.gene_names(i).Data;
end
data.gene_names = gene_names;
clear gene_names 

% create the data_matrix
data_matrix = zeros(data.shape', 'int32');
j = 1;
for i = 1:length(data.data)
    if i == data.indptr(j+1)
       j = j+1;
    end
    data_matrix(data.indices(i), j) = data.data(i);
end
data_matrix = (double(data_matrix'));
data.shape = size(data_matrix);
fprintf('samples = %d\n', data.shape(1));
fprintf('genese = %d\n', data.shape(2));
fprintf('data shape = [%d,%d]\n', data.shape);
fprintf('=============================================================\n');

%% Preprocessing
fprintf('Removing the genes with 0 count from the dataset...\n');
% find genes with zero counts
gene_count = sum(data_matrix);
zero_genes = find(gene_count == 0);
fprintf('%d genes found and removed!\n', length(zero_genes));
% remove genes with zero counts
data_matrix(:, zero_genes) = [];
data.gene_names(zero_genes) = [];
data.genes(zero_genes) = [];
data.shape = size(data_matrix);
fprintf('data shape = [%d,%d]\n', data.shape);
fprintf('=============================================================\n');
% remove genes below cut-off threshold
if  exist('cutoff_thres', 'var')
    fprintf('Low-expression gene filtering...\n');
    fprintf('Removing genes with sum expressed count < %d.\n', cutoff_thres);
    % set cut-off threshold to filter low expression genes
    gene_count = sum(data_matrix);
    cutoff_genes = find(gene_count < cutoff_thres);
    % remove genes below threshold
    data_matrix(:, cutoff_genes) = [];
    data.gene_names(cutoff_genes) = [];
    data.genes(cutoff_genes) = [];
    data.shape = size(data_matrix);
    fprintf('data shape = [%d,%d]\n', data.shape);
end
fprintf('=============================================================\n');
data.data_matrix = data_matrix;
%% Dimension Reduction:
switch dim_red_method
    case 'PCA'
        fprintf('Reducing the dimension of the data using PCA...\n');
        [U, S] = pca(data_matrix);
        data_matrix = projectData(data_matrix, U, red_dim);
    case 'SNE'
        fprintf('Reducing the dimension of the data using SNE...\n');
        data_matrix = tsne(data_matrix, [], red_dim);
    otherwise
        error('Unexpected dimension reduction type! Use PCA or SNE')
end

data.shape = size(data_matrix);
fprintf('data shape = [%d,%d]\n', data.shape);
fprintf('=============================================================\n');

% figure,
% scatter3(data_matrix(:,1), data_matrix(:,2), data_matrix(:,3), 'filled');
% title('Visualization of samples in 3D')

%% Normalize [0 - 1]
if normalize
    fprintf('Normalizing the features to [0-1]...\n');
    data_matrix = bsxfun(@rdivide, ...
                         bsxfun(@minus, data_matrix, min(data_matrix)),...
                         max(data_matrix) - min(data_matrix)...
                         );
    fprintf('=============================================================\n');
end

%% DPMM:
fprintf('Running Non-Parametric clustering...\n');
N = size(data_matrix, 1);
T = 20;
Phi = 1/T*ones(N,T);
alpha = rand*10;
% Run EM!
removed_samples = [];
for i = 1:100
    [gamma,mu_0,lambda,W,nu] = Mstep(data_matrix,Phi,alpha);
    [Phi,alpha] = Estep(data_matrix,gamma,mu_0,lambda,W,nu);
    t = any(isnan(Phi),2);    % find the samples with NaN
    Phi = Phi(~t, :);         % remove element from phi
    data_matrix = data_matrix(~t, :); % remove element from the samples
    data.barcodes = data.barcodes(~t);
    data.shape = size(data_matrix);
    fprintf('iter: %d\tdata shape: [%d,%d]\n', i, data.shape);
end

% find the highest prob of each cluster. i.e. cluster label
[~, Z_hat] = max(Phi,[],2);
% sort the labels from 1 to last found cluster
labels = zeros(1,length(Z_hat));
UZ = sort(unique(Z_hat));
for i = 1:length(UZ)
    labels(Z_hat==UZ(i)) = i;
end
data.labels = labels';
num_labels = length(unique(labels));    % number of labels (derived clusters from DPMM)
fprintf('Clustering finished with %d clusters\n', num_labels);

% find counts of samples in each cluster
fprintf('cluster \t count\n');
fprintf('======= \t =====\n');
cluster_count = zeros(num_labels,1);
for i = 1:num_labels
    cluster_count(i) = sum(labels == i);
    fprintf('%d \t\t\t %d\n', i, cluster_count(i));
end
fprintf('=============================================================\n');
%% Visualization:
fprintf('Visualizing results using t-SNE...\n');
figure('units','normalized','outerposition',[0 0 1 1]);
cl_map_clust = hsv (num_labels);                                            % create a color map based on different clusters
tsne_labels = cl_map_clust(labels,:);
if data.shape(2) ~= 3
    tsne_feats = tsne(data_matrix, tsne_labels, 3);
else
    tsne_feats = data_matrix;
    scatter3(tsne_feats(:,1), tsne_feats(:,2), tsne_feats(:,3), ...
             [], tsne_labels, 'filled')
end
colormap(cl_map_clust); axis on;
ticks=cell(num_labels,1);
for i = 1:num_labels
    ticks{i} = sprintf('%d - %d cells', i, cluster_count(i));
end
colorbar('Ticks', linspace(1/(num_labels * 2),...
                           1 - 1/(num_labels * 2),...
                           num_labels), ...
         'TickLabels', ticks);
title('Visualizing t-SNE results')

uicontrol('style', 'text',...
          'Units', 'normalized',...
          'Position', [.05 .915 .07 .02],...
          'String','Genes Names (Space Separate)');
uicontrol('style', 'edit',...
          'Tag', 'gene_name',...
          'Units', 'normalized',...
          'Position', [.05 .9 .05 .02]);

uicontrol('style', 'pushbutton',...
          'Units', 'normalized',...
          'Position', [.1 .9 .02 .02],...
          'String','go',...
          'Callback',{@go_Callback, data});

fprintf('=============================================================\n');
%%
fprintf('pipeline finished.\n');
toc