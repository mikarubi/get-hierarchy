function aux = auxiliaries

aux.get_hierarchies           = @get_hierarchies;
aux.process_fly               = @process_fly;
aux.process_mouse             = @process_mouse;
aux.add_lattice               = @add_lattice;
aux.core_indices              = @core_indices;

end

function [MC,MR,Hier,MCh,Knot] = get_hierarchies(W,gamm_range,Ih1,Ih2)

if nargin>1
    modu_iter = 1000;
    modu_thr  = 4;
    
    [MC,MR] = multiscale_community(W, modu_iter, gamm_range, Ih1, Ih2);
else
    MC = W;
end

if nargout>1
    MCh = MC;
    Knot = false(1,size(MC,2));
    while 1
        [Hier, MCh, Knot_new] = clean_hierarchies(MCh,modu_thr);
        
        if any(Knot_new)
            Knot = Knot | Knot_new;
            MCh  = MC;
            MCh(:,Knot) = 0;
        else
            break
        end
    end
end

end

function [W,B,Dtwm,Deuc,Ih1,Ih2,ROI_nm,ROI_fullname] = process_fly

load DataSets/flybrain_network_construction Wabs IdxHemi1 IdxHemi2
W = Wabs;
n = size(W,1);
W(isnan(W))  = 0;
W(1:n+1:end) = 0;

B = W;
Ih1 = false(1,n);
Ih2 = false(1,n);
Ih1(IdxHemi1) = 1;
Ih2(IdxHemi2) = 1;

load DataSets/flybrain_tract_lengths VD_medi1 ED_medi0 ROI_nm X Y Z

idx_d1  = find(VD_medi1>0);
mi      = ED_medi0(idx_d1);
di      = VD_medi1(idx_d1);
krnl    = exp((-(bsxfun(@minus,mi,mi.')).^2)/20);
ph0     = -inf;
ph      =  inf;
wt      =    1;
while abs(ph-ph0)>1e-4;
    ph0  = ph;
    ph   = (sqrt(wt).*mi) \ (sqrt(wt).*di);
    se2  = (di-mi*ph).^2;
    se2s = sum(bsxfun(@times,krnl,se2))./sum(krnl);
    wt   = 1./se2s.';
end
Dtwm       = VD_medi1;
idx_d0  = find(~(Dtwm>0) & ~eye(n));
Dtwm(idx_d0) = ph*ED_medi0(idx_d0);

Dtwm = Dtwm*(4*0.6221)/1000;      % 0.6221 is downsampled voxel size

% get centroids and compute Euclidean distance
Cx = cellfun(@median, X);
Cy = cellfun(@median, Y);
Cz = cellfun(@median, Z);

Deuc = sqrt( ...
    bsxfun(@minus, Cx, Cx.').^2 +   ...
    bsxfun(@minus, Cy, Cy.').^2 +   ...
    bsxfun(@minus, Cz, Cz.').^2 );
Deuc = Deuc*(4*0.6221)/1000;      % 0.6221 is downsampled voxel size

nm_map = { ...
    'AL',   'Antennal Lobe'
    'AMMC', 'Antennal Mechanosensory and Motor Center'
    'CCP',  'Caudalcentral Protocerebrum'
    'CMP',  'Caudalmedial Protocerebrum'
    'CVLP', 'Caudal Ventrolateral Protocerebrum'
    'DLP',  'Dorsolateral Protocerebrum'
    'DMP',  'Dorsomedial Protocerebrum'
    'EB',   'Ellipsoid Body'
    'FB',   'Fanshaped Body'
    'IDFP', 'Inferior Dorsofrontal Protocerebrum'
    'LH',   'Lateral Horn'
    'LOB',  'Lobulla'
    'LOP',  'Lobulla Plate'
    'MB',   'Mushroom Body'
    'MED',  'Medulla'
    'PAN',  'Proximal Antennal Protocerebrum'
    'PCB',  'Protocerebral Bridge'
    'SDFP', 'Superior Dorsofrontal Protocerebrum'
    'SOG',  'Subesophageal Ganglion'
    'SPP',  'Superpenduncular Protocerebrum'
    'VLP-D','Ventrolateral Protocerebrum, Dorsal part'
    'VLP-V','Ventrolateral Protocerebrum, Ventral part'
    'VMP',  'Ventromedial Protocerebrum'
    'NOD',  'Noduli'
    'OG',   'Optic Glomerulus'
    'OPTU', 'Optic Tubercle'
    };

ROI_fullname = cellfun(@(i) nm_map{strcmpi(nm_map(:,1),i),2}, ROI_nm, 'uniformoutput',false);

end

function [W,B,Dtwm,Deuc,Ih1,Ih2,ROI_nm,ROI_fullname] = process_mouse

mouse_folder = 'C:/Users/mrubi/Documents/Research/past_projects/2013_mouse/';
addpath(mouse_folder)
load([mouse_folder                  'Data_BMU_20140910.mat'],'Prcl_ID','Prcl_XYZ','dim3','S')
load([mouse_folder 'Network_Analysis_BMU_20140910_-2_4.mat'],'MC_emp','w_emp','b_emp','prcl_idx','D')

W = w_emp;
Dtwm = D{3}/10;                % convert to mm
B = b_emp;

n = size(W,1);
h = n/2;
Ih1 = false(1,n);
Ih2 = false(1,n);
Ih1(   1:h)  = 1;
Ih2(h+(1:h)) = 1;

% get centroids and compute Euclidean distance
Cx = zeros(n,1);
Cy = zeros(n,1);
Cz = zeros(n,1);
for i=1:n
    [cx,cy,cz] = ind2sub(dim3, Prcl_XYZ{prcl_idx(i)});
    Cx(i) = median(cx);
    Cy(i) = median(cy);
    Cz(i) = median(cz);
end

Deuc = sqrt( ...
    bsxfun(@minus, Cx, Cx.').^2 +   ...
    bsxfun(@minus, Cy, Cy.').^2 +   ...
    bsxfun(@minus, Cz, Cz.').^2 );
Deuc = Deuc/10;                % convert to mm

ROI_nm       = arrayfun( @(i) S{4}{S{1}==i}, Prcl_ID(prcl_idx), 'uniformoutput', false);
ROI_fullname = arrayfun( @(i) S{3}{S{1}==i}, Prcl_ID(prcl_idx), 'uniformoutput', false);

end

function [NULL1,BAND1] = add_lattice(NULL1,BAND1,W,B,D,Ih1,Ih2)

n = size(W,1);
eidx = find(~eye(n));

j = size(NULL1,2)+1;
medd = median(D(eidx));
for i=1:size(NULL1,1);
    R = rand(n);
    R(Ih2,Ih2) = R(Ih1,Ih1);
    R(Ih2,Ih1) = R(Ih1,Ih2);
    
    [~,ord_d] = sort(D(eidx)+R(eidx)*medd/100);
    [~,ord_w] = sort(W(eidx),'descend');
    
    WL = zeros(n);
    WL(eidx(ord_d)) = W(eidx(ord_w));
    
    BL = zeros(n);
    BL(eidx(ord_d)) = B(eidx(ord_w));
    
    NULL1{i,j} = WL;
    BAND1{i,j} = BL;
end
end

function Core = core_indices(W, Ih1, Ih2)

Jh1 = find(Ih1);
Jh2 = find(Ih2);

n = length(W);
Wh = W;
Wh(Ih1,:) = Wh(Ih1,:)+Wh(Ih2,:);
Wh(:,Ih1) = Wh(:,Ih1)+Wh(:,Ih2);
Wh = Wh(Ih1,Ih1);
[~, Ord]  = sort(sum(Wh)+sum(Wh,2).');

Core = false(1, n);
switch n;
    case 112; Core([Jh1(Ord(50:end)) Jh2(Ord(50:end))]) = 1;
    case 49;  Core([Jh1(Ord(21:end)) Jh2(Ord(21:end))]) = 1;
end

end
