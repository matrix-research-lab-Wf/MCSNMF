function run_all_proposed_method()
%RUN_ALL_PROPOSED_METHOD Reproduce MCS-NMF on all fixed paired splits.
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
datasets={'YaleB','COIL20','Optdigits','PIE','AR','MSRA25','COIL100','MNIST','Palm'};
for i = 1:numel(datasets)
    dataset = datasets{i};
    split_dir = fullfile(root,'protocol',dataset);
    splits = dir(fullfile(split_dir,[dataset '_seed*.mat']));
    for j = 1:numel(splits)
        result = run_mcsnmf(fullfile(root,'data',[dataset '.mat']), ...
            fullfile(split_dir,splits(j).name),2,3);
        out_dir = fullfile(root,'results','reproduction',dataset);
        if ~exist(out_dir,'dir'), mkdir(out_dir); end
        save(fullfile(out_dir,[splits(j).name(1:end-4) '_MCSNMF.mat']),'result');
    end
end
end
