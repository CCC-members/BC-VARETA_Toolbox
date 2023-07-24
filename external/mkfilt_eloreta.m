function [A,Wout]=mkfilt_eloreta(L,regu,W)
% makes spatial filter according to eLoreta
% usage  A=mkfilt_eloreta(L); or  A=mkfilt_eloreta(L,regu);
%
% input L:  NxMxP leadfield tensor for N channels, M voxels, and
%           P dipole directions. Typically P=3. (If you do MEG for
%           a spherical volume conductor or reduce the rank, you must
%           reduce L such that it has full rank for each voxel, such that,
%           e.g., P=2)
%       regu: optional regularization parameter (default is .05 corresponding
%             to 5% of the average of the eigenvalues of some matrix to be inverted.)
%
% output A: NxMxP tensor of spatial filters. If x is the Nx1 data vector at time t.
%           then A(:,m,p)'*x is the source activity at time t in voxel m in source direction
%           p.
%
% code implemented in the same style as mkfilt_eloreta by Guido Nolte
% please cite
% “R.D. Pascual-Marqui: Discrete, 3D distributed, linear imaging methods of electric neuronal activity. Part 1: exact, zero
% error localization. arXiv:0710.3341 [math-ph], 2007-October-17, http://arxiv.org/pdf/0710.3341 ”

if nargin<2
    regu =.05;
end

[nchan, ng, ndum] = size(L);
LL=zeros(nchan,ndum,ng);
for i = 1:ndum
    LL(:,i,:) = L(:,:,i);
end
LL = reshape(LL,nchan,ndum*ng);

u0 = eye(nchan);
if nargin<3
    W = reshape(repmat(eye(ndum),1,ng),ndum,ndum,ng);
end

Winv = zeros(ndum,ndum,ng);
winvkt = zeros(ng*ndum,nchan);
for i = 1:ng
    Winv(:,:,i) = (inv(W(:,:,i) + trace(W(:,:,i))/(ndum*10^6)));
end
for i = 1:ng
    winvkt(ndum*(i-1)+1:ndum*i,:) = Winv(:,:,i)*LL(:,ndum*(i-1)+1:ndum*i)';
end
kwinvkt = LL*winvkt;
M =inv(kwinvkt + regu*u0);

for i = 1:ng
    Lloc = squeeze(L(:,i,:));
    Mb = Lloc'*M*Lloc;
    Mb = Mb+eye(ndum)*trace(Mb)/10^10;
    W(:,:,i) = sqrtm(Mb);
end

ktm=LL'*M;
A=zeros(nchan,ng,ndum);

for i=1:ng;
    A(:,i,:)=(Winv(:,:,i)*ktm(ndum*(i-1)+1:ndum*i,:))';
end
Wout=W;
return
