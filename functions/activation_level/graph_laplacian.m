function [D,D3D] = graph_laplacian(faces,s)
Nv  = max(faces(:));
D   = zeros(Nv,Nv);
D3D = zeros(3*Nv,3*Nv);
for node = 1:Nv
    [row,col]               = find(faces == node);
    index                   = faces(row,:);
    index                   = index(:);
    index                   = setdiff(index,node);
    D(node,node)            = size(index,1);
    D(node,index)           = -1;
    D(index,node)           = -1;
    D3D(3*node,3*node)      = size(index,1);
    D3D(3*node-1,3*node-1)  = size(index,1);
    D3D(3*node-2,3*node-2)  = size(index,1);
    D3D(3*node,3*index)     = -1;
    D3D(3*node-1,3*index-1) = -1;
    D3D(3*node-2,3*index-2) = -1;
    D3D(3*index,3*node)     = -1;
    D3D(3*index-1,3*node-1) = -1;
    D3D(3*index-2,3*node-2) = -1;
end
D    = sparse(D);
Deg  = spdiags(diag(D),0,Nv,Nv);
A    = abs(D - Deg);
Ia   = speye(length(A));
D    = Ia - s*A + (s^2)*(Deg - Ia);
D3D  = sparse(D3D);
Deg  = spdiags(diag(D3D),0,3*Nv,3*Nv);
A    = abs(D3D - Deg);
Ia   = speye(length(A));
D3D  = Ia - s*A + (s^2)*(Deg - Ia);
end