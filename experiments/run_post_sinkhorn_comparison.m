function run_post_sinkhorn_comparison(dataset)
%RUN_POST_SINKHORN_COMPARISON Paired Row-only/Post-Sinkhorn/MCS-NMF runs.
here=fileparts(mfilename('fullpath')); root=fileparts(here); addpath(fullfile(root,'code'));
[p,k]=parameters(dataset); [X,y]=load_data(fullfile(root,'data',[dataset '.mat']));
outdir=fullfile(root,'results','post_sinkhorn_reproduction');if~exist(outdir,'dir'),mkdir(outdir);end
outfile=fullfile(outdir,[dataset '.csv']);
if exist(outfile,'file')~=2
 fid=fopen(outfile,'w');fprintf(fid,'dataset,seed,variant,ACC,NMI,Purity,Fscore,Precision,Recall,iterations,runtime_seconds\n');fclose(fid);
end
T=readtable(outfile);done=[];if~isempty(T),done=unique(T.seed(strcmp(T.variant,'MCS-NMF')));end
for seed=20260627:20260646
 if any(done==seed),continue;end
 S=load(fullfile(root,'protocol',dataset,sprintf('%s_seed%d.mat',dataset,seed)));
 lab=logical(S.lab(:));un=~lab;order=[find(lab);find(un)];Xo=X(:,order);yo=y(order);
 n=numel(y);c=max(y);nl=sum(lab);Y=zeros(n,c);Y(sub2ind([n c],(1:nl)',yo(1:nl)))=1;
 L=build_high_order_laplacian(Xo,p,k);labo=[true(nl,1);false(n-nl,1)];
 base=struct('alpha',1000,'maxit',100,'tol',1e-7,'seed',seed);
 rowopt=base;rowopt.mass_conservation=false;t=tic;[~,Frow,irow]=mcsnmf(Xo,L,Y,labo,rowopt);rtrow=toc(t);
 target=(n-nl)*sum(Y(1:nl,:),1)/nl;Fpost=Frow;Fpost(nl+1:end,:)=transport_project(Frow(nl+1:end,:),target,300);
 fullopt=base;fullopt.mass_conservation=true;t=tic;[~,Ffull,ifull]=mcsnmf(Xo,L,Y,labo,fullopt);rtfull=toc(t);
 variants={'Row-only','Post-Sinkhorn','MCS-NMF'};Fs={Frow,Fpost,Ffull};infos={irow,irow,ifull};rts=[rtrow 0 rtfull];
 fid=fopen(outfile,'a');
 for v=1:3
  [~,po]=max(Fs{v},[],2);pred=zeros(n,1);pred(order)=po;m=metrics6(y(un),pred(un));
  fprintf(fid,'%s,%d,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.9g\n',dataset,seed,variants{v},m.ACC,m.NMI,m.Purity,m.Fscore,m.Precision,m.Recall,infos{v}.iterations,rts(v));
 end
 fclose(fid);fprintf('%s seed=%d complete\n',dataset,seed);
end
end
function [p,k]=parameters(name)
switch name
 case {'YaleB','COIL20','COIL100'},p=2;k=3;
 case 'PIE',p=3;k=5;
 case {'AR','Optdigits','MNIST'},p=5;k=1;
 case {'MSRA25','Palm'},p=10;k=1;
 otherwise,error('Unknown dataset %s',name);end
end
function [X,y]=load_data(path)
D=load(path);if isfield(D,'X'),X=double(D.X);else,X=double(D.fea)';end
if isfield(D,'y'),y=D.y(:);else,y=D.gnd(:);end;if size(X,2)~=numel(y),X=X';end
[~,~,y]=unique(y);y=double(y);X=X+max(0,-min(X(:)));X=bsxfun(@rdivide,X,max(sqrt(sum(X.^2,1)),eps));
end
function P=transport_project(K,target,maxit)
P=max(K,1e-300);target=max(target,1e-12);target=target*(size(P,1)/sum(target));
for it=1:maxit,P=bsxfun(@rdivide,P,max(sum(P,2),1e-300));P=bsxfun(@times,P,target./max(sum(P,1),1e-300));end
end
function m=metrics6(y,p)
[~,~,y]=unique(y(:));[~,~,p]=unique(p(:));C=accumarray([y p],1,[max(y) max(p)]);n=sum(C(:));m.Purity=sum(max(C,[],1))/n;
tp=sum(sum(C.*(C-1)))/2;pp=sum(sum(C,1).*(sum(C,1)-1))/2;tt=sum(sum(C,2).*(sum(C,2)-1))/2;
m.Precision=tp/max(pp,eps);m.Recall=tp/max(tt,eps);m.Fscore=2*m.Precision*m.Recall/max(m.Precision+m.Recall,eps);
pa=sum(C,2)/n;pb=sum(C,1)/n;q=C/n;I=0;for i=1:size(C,1),for j=1:size(C,2),if q(i,j)>0,I=I+q(i,j)*log(q(i,j)/(pa(i)*pb(j)));end,end,end
Ha=-sum(pa(pa>0).*log(pa(pa>0)));Hb=-sum(pb(pb>0).*log(pb(pb>0)));m.NMI=I/max(sqrt(Ha*Hb),eps);m.ACC=mean(bestMap(y,p)==y);
end
