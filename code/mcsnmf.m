function [U,F,info] = mcsnmf(X,L,Y,labeled,opts)
%MCSNMF Mass-conserving semantic nonnegative matrix factorization.
% X is d-by-n, L is n-by-n, Y is n-by-c, and labeled is n-by-1.
if nargin < 5, opts = struct(); end
alpha = get_opt(opts,'alpha',1000);
maxit = get_opt(opts,'maxit',100);
tol = get_opt(opts,'tol',1e-7);
seed = get_opt(opts,'seed',1);
rng(seed,'twister');
[d,n] = size(X); c = size(Y,2);
unlabeled = ~logical(labeled(:));
target = sum(unlabeled) * sum(Y(labeled,:),1) / sum(labeled);
U = max(rand(d,c),eps);
F = max(rand(n,c),eps);
F(labeled,:) = Y(labeled,:);
F(unlabeled,:) = transport_project(F(unlabeled,:),target,300);
objective = zeros(maxit+1,1);
objective(1) = obj(X,U,F,L,alpha);
step = [];
for it = 1:maxit
    U = U .* ((X*F) ./ max(U*(F'*F),1e-12));
    grad = 2*(F*(U'*U)-X'*U+alpha*(L*F));
    base = obj(X,U,F,L,alpha);
    if isempty(step), step = 1/max(2*(norm(U'*U,2)+alpha),1e-12);
    else, step = min(1,1.5*step); end
    accepted = false;
    for bt = 1:40
        K = max(F(unlabeled,:),1e-12).*exp(max(min(-step*grad(unlabeled,:),50),-50));
        Ft = F; Ft(unlabeled,:) = transport_project(K,target,200); Ft(labeled,:) = Y(labeled,:);
        if obj(X,U,Ft,L,alpha) <= base, F = Ft; accepted = true; break; end
        step = step/2;
    end
    if ~accepted, step = 0; end
    objective(it+1) = obj(X,U,F,L,alpha);
    if it >= 30 && abs(objective(it+1)-objective(it))/max(1,abs(objective(it))) <= tol, break; end
end
objective = objective(1:it+1);
Fu = F(unlabeled,:);
info.iterations = it;
info.objective = objective;
info.max_row_error = max(abs(sum(Fu,2)-1));
info.max_column_error = max(abs(sum(Fu,1)-target));
end

function P = transport_project(K,target,maxit)
P = max(K,1e-300); target = max(target,1e-12);
target = target*(size(P,1)/sum(target));
for t = 1:maxit
    P = bsxfun(@rdivide,P,max(sum(P,2),1e-300));
    P = bsxfun(@times,P,target./max(sum(P,1),1e-300));
    if mod(t,10)==0 && max(abs(sum(P,2)-1))<1e-10 && max(abs(sum(P,1)-target))<1e-8, break; end
end
end

function value = obj(X,U,F,L,alpha)
value = sum(X(:).^2)-2*sum(sum((U'*X).*F')) ...
    + sum(sum((U'*U).*(F'*F))) + alpha*sum(sum(F.*(L*F)));
end

function value = get_opt(opts,name,default)
if isfield(opts,name), value = opts.(name); else, value = default; end
end
