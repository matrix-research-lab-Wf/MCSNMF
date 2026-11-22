function run_all_proposed_method()
%RUN_ALL_PROPOSED_METHOD Reproduce all proposed-method and mechanism runs.
datasets={'YaleB','COIL20','Optdigits','PIE','AR','MSRA25','COIL100','MNIST','Palm'};
for i=1:numel(datasets),run_post_sinkhorn_comparison(datasets{i});end
end
