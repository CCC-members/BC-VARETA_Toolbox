function [A,Wout]=mkfilt_mne(L,alpha,C)
% [A]=mkfilt_lcmv(L,C,alpha)
%
% Input:
% L: Lead field matrix, NxMx3 matrix for N channels, M voxels and 3 dipole
%    directions
% C:  NxN matrix for N channels covariance matrix or real part of cross-spectrum
%    (The program uses only the real of C)
% 
% Output
% A  : 3-dimension filter
% 
% History
% code implemented in the same style than mkfilt_lcmv and mkfilt_eloreta by Guido Nolte
% please cite as 
% Hämäläinen, M.S., Ilmoniemi, R.J. Interpreting magnetic fields of the brain: minimum norm estimates. Med. Biol. Eng. Comput. 32, 35–42 (1994). https://doi.org/10.1007/BF02512476

% License
%   Copyright (C) 2015
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see http://www.gnu.org/licenses/.

[nchan,ns,ndim] = size(L);
L = reshape(L,nchan,ns*ndim);

if nargin<2
    alpha=.05*trace(C)/length(C);
end

if nargin<3
    C = eye(nchan);
end
C = real(C);
[U,D] = eig(C);
C = diag(abs(diag(D)).^(-1/2))*U';
L = C*L;

A = (L'/(L*L' + alpha))*C;
Wout = alpha*ones(ns,1);