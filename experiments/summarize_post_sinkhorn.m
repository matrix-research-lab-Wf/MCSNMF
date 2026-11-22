function summarize_post_sinkhorn()
%SUMMARIZE_POST_SINKHORN Aggregate paired runs with percentile-bootstrap CIs.
here=fileparts(mfilename('fullpath'));root=fileparts(here);indir=fullfile(root,'results','post_sinkhorn_reproduction');
files=dir(fullfile(indir,'*.csv'));rows={};rng(20261123,'twister');
for q=1:numel(files)
 T=readtable(fullfile(files(q).folder,files(q).name));ds=files(q).name(1:end-4);
 seeds=unique(T.seed);da=zeros(numel(seeds),1);dn=da;
 for i=1:numel(seeds)
  m=T(T.seed==seeds(i)&strcmp(T.variant,'MCS-NMF'),:);p=T(T.seed==seeds(i)&strcmp(T.variant,'Post-Sinkhorn'),:);
  assert(height(m)==1&&height(p)==1);da(i)=m.ACC-p.ACC;dn(i)=m.NMI-p.NMI;
 end
 [la,ua]=bootci(da,10000);[ln,un]=bootci(dn,10000);
 rows(end+1,:)={ds,numel(seeds),mean(da),la,ua,sum(da>0),sum(da<0),sum(da==0),mean(dn),ln,un,sum(dn>0),sum(dn<0),sum(dn==0)}; %#ok<AGROW>
end
R=cell2table(rows,'VariableNames',{'dataset','n','MCS_minus_Post_ACC','ACC_CI_low','ACC_CI_high','ACC_wins','ACC_losses','ACC_ties','MCS_minus_Post_NMI','NMI_CI_low','NMI_CI_high','NMI_wins','NMI_losses','NMI_ties'});
writetable(R,fullfile(root,'results','post_sinkhorn_reproduced_summary.csv'));disp(R);
end
function [lo,hi]=bootci(x,B)
n=numel(x);v=zeros(B,1);for b=1:B,v(b)=mean(x(randi(n,n,1)));end;v=sort(v);lo=v(max(1,round(.025*B)));hi=v(min(B,round(.975*B)));
end
