function [Hier, MCh, Knot] = clean_hierarchies(MC,thr)

[~, n] = size(MC);                                          % number of hierarchies and nodes
[~, MI] = partition_distance_matrix(MC);

MI(abs(MI-1) >1e-5) = 0;
MI(abs(MI-1)<=1e-5) = 1;                                    % binarize distance
if ~isequalwithequalnans(MI,MI.')
    error('MI symmetry error.')
end

[~, I, U] = unique(-MI,'rows');                             % unique hiearchies
IdxU = cell(size(I));
for i = 1:numel(I)
    IdxU{i} = find(U==U(I(i)));
end

Hier = IdxU(cellfun(@numel,IdxU)>=thr);                     % hierarchical indices
MCh  = MC(cellfun(@(i) i(1), Hier), :);                     % hierarchical partitions

Knot = false(1, size(MC,2));
for h = 2:numel(Hier)
    Uh = nonzeros(unique(MCh(h,:)));                        % unique modules in hierarchy
    
    for i = 1:numel(Uh)                                     % loop over unique modules
        idx_ui = MCh(h,:)==Uh(i);
        Uh_prev_ui = unique(MCh(h-1,idx_ui));               % get module in prev hierarchy
        
        if numel(Uh_prev_ui)>1
            size_Uh_prev_ui = arrayfun(@(j) nnz(MCh(h-1,idx_ui) == j), Uh_prev_ui);
            
            [~,idxmax_Uh_prev_ui] = max(size_Uh_prev_ui);   % get size of modules
            for j = [1:(idxmax_Uh_prev_ui-1) (idxmax_Uh_prev_ui+1):numel(Uh_prev_ui)]
                idx_knot = idx_ui & (MCh(h-1,:)==Uh_prev_ui(j));
                MCh(:, idx_knot) = 0;                       % remove in MCh
                Knot(  idx_knot) = 1;                       % assign to hi-par
            end
        end
        
        idx_ui = MCh(h,:)==Uh(i);
        if (nnz(idx_ui)/n)<0.02                             % put small mouse modules in knot
            MCh(:,idx_ui)=0;                                % remove in MCh
            Knot( idx_ui)=1;                                % assign to hi-par
        end
    end
end
