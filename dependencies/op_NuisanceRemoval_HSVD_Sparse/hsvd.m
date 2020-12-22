
function [z,f,Z,a]= hsvd(xt, L, N, n, dt)
% n: number of most dominant sing. values to keep; MUST <=N, model
% order.

% Compute data matrix
c=xt(1:(L+1),1).';
r=xt((L+1):N,1).';
X=hankel(c,r);

% take SVD 
[U,S,V] = svd(X);

% discarding noise subspace: leave only n_singval principal sing. values/vectors
U=U(:,1:n);
V=V(1:n,:);

% finding U1=U_ - omitting the last row of U containing N largest sing.
% vectors
g=size(U,1);
U1=U(1:g-1,:);
 
% finding U2=U- - omitting the first row of U containing N largest sing.
% vectors
U2=U(2:g,:); 

z=eig(pinv(U1)*U2);
% or
% u=U(size(U,1),:);
% Z=(eye(n)+(u'*u)/(1-u*u'))*U1'*U2;
% z=eig(Z);

% Sorting z by angle(z)
[~,order]=sort(angle(z));                      
z=z(order);
f = angle(z)/(dt*2*pi);

% finding a's through already sorted z
% A=Matrix_z(z.',N,M);
% [U,S,V]=svd([A xt]);     % S: singular values in decreasing order
% a=-(1/V(N+1,N+1))*V(:,N+1);          % take right sing. vector corresponding to smallest sing. value and normalization
% a=a(1:N,1);           % throw last term, it corresponds to z^n in the
%a=pinv(Matrix_Z_new(z.',N,M))*xt;                 % found a's must be already in the correct order compared with original already sorted a's

Z=ones(length(xt),length(z));
temp=reshape(z, 1, length(z));
for i=2:length(xt) %generate Vandermonde matrix B based on poles
    Z(i,:)=Z(i-1,:).*temp;
end

if nargout>3
    a=Z\xt;
end


