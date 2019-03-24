function [VIn,MIn] = partition_distance_matrix(CX,CY)

%#ok<*ASGLU>

s = (nargin==1);
if s
    CY = CX;
end

d = 10.^ceil(log10(double(max(max([CX; CY]))+1)));

n = size(CX,2);                                                 %n
x = size(CX,1);
CX = uint32(CX);
HX = zeros(x,1);
for i = 1:x
    Px = nonzeros(histif(CX(i,:)))/n;                           %P(x)
    HX(i) = -sum(Px.*log(Px));                                  %H(x)
end

assert(n == size(CY,2));
y = size(CY,1);
CY = uint32(CY);
HY = zeros(y,1);
for j = 1:y
    Py = nonzeros(histif(CY(j,:)))/n;                           %P(y)
    HY(j) = -sum(Py.*log(Py));                                  %H(y)
end

VIn = zeros(x,y);
MIn = zeros(x,y);
for i = 1:x
    j_idx = (s*(i-1)+1):y;
    for j = j_idx
        Pxy = nonzeros(histif(d*CX(i,:) + CY(j,:)))/n;          %P(x,y)
        Hxy = -sum(Pxy.*log(Pxy));                              %H(x,y)
        VIn(i, j) = (2*Hxy - HX(i) - HY(j))/log(n);             %VIn
        MIn(i, j) = 2*(HX(i) + HY(j) - Hxy)/(HX(i) + HY(j));    %MIn
    end
    if s
        VIn(j_idx, i) = VIn(i, j_idx);
        MIn(j_idx, i) = MIn(i, j_idx);
    end
end

