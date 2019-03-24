function [MC, MR] = multiscale_community(W, iter, gamm_range, Ih1, Ih2)

% MC: vector of community affiliations across gammas
% MR: raw (unprocessed) vector of community affiliations across gammas

if iscell(W)
    l = numel(W);
    n = size(W{1}, 2);
else
    n = size(W, 2);
end

MC = zeros(numel(gamm_range),n);
MR = zeros(numel(gamm_range),n);
for i=1:numel(gamm_range);
    gamm = gamm_range(i);
    
    if iscell(W)
        w0 = zeros(n);
        for j = 1:l
            m0 = community_louvain_mex(W{j}, gamm);
            w0 = w0 + bsxfun(@eq, m0, m0.')/l;
        end
    else
        w0 = W;
    end
    w1 = ~eye(n);
    for h = 1:100
        w0  = w0.*w1;
        w1  = zeros(n);
        MCi = zeros(iter,n);
        MRi = zeros(iter,n);
        for j = 1:iter
            q0 = -inf;
            q1 = 0;
            m1 = 1:n;
            while q1-q0 > 1e-5;
                q0 = q1;
                [m1, q1] = community_louvain_mex(w0, gamm, m1);
            end
            m2 = community_symmetry(m1, Ih1, Ih2);
            
            MRi(j,:) = m1;            
            MCi(j,:) = m2;
            w1 = w1 + bsxfun(@eq, m2, m2.')/iter;
        end
        
        converged = all(abs(w1(:))<1e-10 | abs(w1(:)-1)<1e-10);
        if converged
            break
        end
    end
    
    if converged
        j = 1;
    else
        ctic = tic;
        [~, MI_] = partition_distance_matrix(MCi);
        [~,   j] = max(sum(MI_ > 1-1e-10));

        disp(['enforced convergence (gamm, toc): ' num2str([gamm toc(ctic)]) '.'])
    end
    
    MR(i,:) = MRi(j,:);
    
    MCi(j,Ih2) = MCi(j,Ih1);
    [~, ~, MC(i,:)] = unique(MCi(j,:));    
end

end

function Msym = community_symmetry(M, Ih1, Ih2)

Msym = M(:).';

Ma  = Msym(Ih1);
Mb  = Msym(Ih2);
Jac = zeros(max(Ma),max(Mb));                   % Inter-community Jaccard
for a=1:max(Ma);
    for b=1:max(Mb);
        Jac(a,b) =         ...
            nnz(Ma==a & Mb==b)/   ...
            nnz(Ma==a | Mb==b);
    end
end

uneqMod = Jac>1e-5 & Jac<(1-1e-5);              % Jaccard neither 0 nor 1

if any(any(uneqMod))
    [ma,mb] = find(uneqMod);                    % find all unequal modules
    ma = ma(:);
    mb = mb(:);                                 % enforce column orientation
    [~,idx] = max(Jac(uneqMod));                % get the most equal module
    
    Ma(any(bsxfun(@eq,Ma,ma),1)) = ma(idx);     % absorb all other modules
    Mb(any(bsxfun(@eq,Mb,mb),1)) = mb(idx);     %   into this module
    
    Msym(Ih1) = Ma;
    Msym(Ih2) = Mb;
end
end
