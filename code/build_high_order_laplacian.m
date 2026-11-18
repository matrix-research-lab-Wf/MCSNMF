function L = build_high_order_laplacian(X,p,k)
%BUILD_HIGH_ORDER_LAPLACIAN Construct the symmetric high-order operator.
% X is d-by-n; p is the neighborhood size; k is the diffusion order.
n = size(X,2); D2 = max(0,bsxfun(@plus,sum(X.^2,1)',sum(X.^2,1))-2*(X'*X));
[~,idx] = sort(D2,2,'ascend'); idx = idx(:,2:min(p+1,n));
H = sparse(repmat((1:n)',1,size(idx,2)),idx,1,n,n);
H = H + speye(n); W = speye(n);
invDe = spdiags(1./max(sum(H,2),eps),0,n,n);
S = H'*W*invDe*H; A = (S+S')/2;
deg = max(sum(A,2),eps); Dh = spdiags(1./sqrt(deg),0,n,n);
P = Dh*A*Dh; P = (P+P')/2;
L = speye(n)-P^k; L = (L+L')/2;
end
