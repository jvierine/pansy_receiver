function tibet_watch(input_folder,output_folder,config)
% Successfully completed sources are skipped on subsequent scans and restarts.
if nargin<3, config=struct(); end
names={'RunOnce','PollSeconds','StableSeconds','Format','PointingFolder'};
values={false,10,20,'auto',''};
for j=1:numel(names), if ~isfield(config,names{j}),config.(names{j})=values{j};end,end
assert(exist(input_folder,'dir')==7,'Input folder does not exist');
if exist(output_folder,'dir')~=7, mkdir(output_folder); end
previous=struct('name',{},'bytes',{},'datenum',{});
while true
    files=dir(fullfile(input_folder,'20*.raw'));
    [~,order]=sort({files.name}); files=files(order);
    for j=1:numel(files)
        file=files(j); [~,stem]=fileparts(file.name);
        target=fullfile(output_folder,stem);
        done=fullfile(target,'catalogue.mat');
        if exist(done,'file')==2
            record=load(done,'result'); r=record.result;
            if r.complete && r.source_bytes==file.bytes && r.source_datenum==file.datenum, continue; end
            warning('Source changed or incomplete result: %s; existing results preserved',file.name); continue;
        end
        if exist(target,'dir')==7 || exist([target '.partial'],'dir')==7
            warning('Incomplete output preserved for inspection: %s',target); continue;
        end
        old=find(strcmp({previous.name},file.name),1);
        stable=~isempty(old) && previous(old).bytes==file.bytes && previous(old).datenum==file.datenum;
        if config.RunOnce, stable=true; end
        if ~stable || (now-file.datenum)*86400<config.StableSeconds, continue; end
        options=struct('Format',config.Format);
        if ~isempty(config.PointingFolder)
            options.PointingFile=fullfile(config.PointingFolder,[stem '_pointing.mat']);
        end
        try
            tibet_process_file(fullfile(input_folder,file.name),target,options);
        catch problem
            warning('Not marked complete: %s: %s',file.name,problem.message);
        end
    end
    if config.RunOnce, break; end
    previous=files;
    pause(config.PollSeconds);
end
end
