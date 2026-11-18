function result = run_mcsnmf(data_file,split_file,p,k)
%RUN_MCSNMF Run one saved paired experiment.
D = load(data_file);
if isfield(D,'X'), X = double(D.X); else, X = double(D.fea)'; end
if isfield(D,'y'), y = D.y(:); else, y = D.gnd(:); end
if size(X,2) ~= numel(y), X = X'; end
[~,~,y] = unique(y); y = double(y);
X = X + max(0,-min(X(:)));
X = bsxfun(@rdivide,X,max(sqrt(sum(X.^2,1)),eps));
S = load(split_file); labeled = logical(S.lab(:));
c = max(y); Y = zeros(numel(y),c);
rows = find(labeled); Y(sub2ind(size(Y),rows,y(rows))) = 1;
L = build_high_order_laplacian(X,p,k);
opts = struct('alpha',1000,'maxit',100,'tol',1e-7);
if isfield(S,'seed'), opts.seed = double(S.seed); end
[U,F,info] = mcsnmf(X,L,Y,labeled,opts);
[~,pred] = max(F,[],2);
result.U = U; result.F = F; result.prediction = pred;
result.labeled = labeled; result.info = info;
result.metrics = clustering_metrics(y(~labeled),pred(~labeled));
end

function m = clustering_metrics(y,p)
[~,~,y] = unique(y(:)); [~,~,p] = unique(p(:));
C = accumarray([y p],1,[max(y) max(p)]); n = sum(C(:));
m.Purity = sum(max(C,[],1))/n;
tp = sum(sum(C.*(C-1)))/2;
predpairs = sum(sum(C,1).*(sum(C,1)-1))/2;
truepairs = sum(sum(C,2).*(sum(C,2)-1))/2;
m.Precision = tp/max(predpairs,eps); m.Recall = tp/max(truepairs,eps);
m.Fscore = 2*m.Precision*m.Recall/max(m.Precision+m.Recall,eps);
pa=sum(C,2)/n; pb=sum(C,1)/n; q=C/n; I=0;
for i=1:size(C,1), for j=1:size(C,2), if q(i,j)>0, I=I+q(i,j)*log(q(i,j)/(pa(i)*pb(j))); end, end, end
Ha=-sum(pa(pa>0).*log(pa(pa>0))); Hb=-sum(pb(pb>0).*log(pb(pb>0)));
m.NMI=I/max(sqrt(Ha*Hb),eps);
m.ACC=mean(bestMap(y,p)==y);
end
