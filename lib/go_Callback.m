function [ output_args ] = go_Callback( src, event, data )
%GO_CALLBACK Summary of this function goes here
%   Detailed explanation goes here
num_labels = max(data.labels);

% extract selected genes from string
genes_str = get(findobj(gcf, 'Tag','gene_name'), 'String');
selected_genes = strsplit(genes_str, ' ');

% create a matrix for heatmap visualization
genes_count = zeros(num_labels, length(selected_genes));
for i = 1:length(selected_genes)
    idx = find(strcmp(data.gene_names, selected_genes{i}));
    if isempty(idx)
        msg = sprintf('Unexpected gene name: %s ', selected_genes{i});
        errordlg(msg)
    end
    for j = 1:num_labels
        genes_count(j, i) = sum(data.data_matrix(data.labels == j, idx));
    end
end

% HMobj = HeatMap(genes_count, ...
%                 'RowLabels', (1:num_labels), ...
%                 'ColumnLabels', selected_genes, ...
%                 'ColumnLabelsRotate', 45);

figure,
HMobj = heatmap(genes_count,...
                selected_genes,...
                1:num_labels,...
                '%i',...
                'TickAngle', 45,...
                'ShowAllTicks', true,...
                'Colorbar', true);
title('Sum of expressed genes clusters')



fprintf('hello\n')
end
