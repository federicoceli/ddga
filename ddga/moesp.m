function [ss,ssfun]=moesp(y,u,d)
% MOESP Multivariable Output Error State Space Approach of Subspace Identification
%   [SS,SSMAT] = MOESP(Y,U,D) identifies the observable subspace based on
%   the measured (N by ny) output matrix, Y with N  sampling points and ny
%   variables, and the corresponding N by nu input matrix, U. The third
%   parameter d is the embeded dimension, which should be larger than the
%   order of the system to be identified. 
%   
%   The function returns the scores of the subspace, SS and a function
%   handle SSMAT, for further identification of state space matrices. 
%
%   The system order can determined from the returned score vector, SS so
%   that sum(SS(1:n)) ~ sum(SS).
%
%   Once the system order is determined, the underline dynamic system is
%   identified by calling the returned function handle:
%
%   [A,B,C,D] = SSMAT(n)
%   
%   to represent a state space model:
%
%       x(k+1) = Ax(k) + Bu(k)
%       y(k)   = Cx(k) + Du(k)
%
% See also: n4sid, subid

% Version 1.0 by Yi Cao at Cranfield University on 27th April 2008

% Input and output check
error(nargchk(1,3,nargin));
error(nargoutchk(0,4,nargout));

[ndat,ny]=size(y);
[mdat,nu]=size(u);
if ndat~=mdat
    error('Y and U have different length.')
end

% block Hankel matrix
N=ndat-d+1;
Y = zeros(d*ny,N);
U = zeros(d*nu,N);
sN=sqrt(N);
sy=y'/sN;
su=u'/sN;
for s=1:d
    Y((s-1)*ny+1:s*ny,:)=sy(:,s:s+N-1);
    U((s-1)*nu+1:s*nu,:)=su(:,s:s+N-1);
end

% LQ decomposition
R=triu(qr([U;Y]'))';
R=R(1:d*(ny+nu),:);

% SVD
R22 = R(d*nu+1:end,d*nu+1:end);
[U1,S1]=svd(R22);

% sigular value
ss = diag(S1);
% n=find(cumsum(ss)>0.85*sum(ss),1);

ssfun = @ssmat;

    function [A,B,C,D]=ssmat(n)
        % C and A
        Ok = U1(:,1:n)*diag(sqrt(ss(1:n)));
        C=Ok(1:ny,:);
        A=Ok(1:ny*(d-1),:)\Ok(ny+1:d*ny,:);

        % B and D
        L1 = U1(:,n+1:end)';
        R11 = R(1:d*nu,1:d*nu);
        R21 = R(d*nu+1:end,1:d*nu);
        M1 = L1*R21/R11;
        m = ny*d-n;
        M = zeros(m*d,nu);
        L = zeros(m*d,ny+n);
        for k=1:d
            M((k-1)*m+1:k*m,:)=M1(:,(k-1)*nu+1:k*nu);
            L((k-1)*m+1:k*m,:)=[L1(:,(k-1)*ny+1:k*ny) L1(:,k*ny+1:end)*Ok(1:end-k*ny,:)];
        end
        DB=L\M;
        D=DB(1:ny,:);
        B=DB(ny+1:end,:);
    end
end


