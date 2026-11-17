function stats=tibet_watch(input_folder,output_folder,config)
% Folder worker: complete files only; date-organized cuts; opt-in verified cleanup.
if nargin<3, config=struct(); end
names={'RunOnce','PollSeconds','StableSeconds','Format','PointingFolder', ...
    'MakeOverview','DeleteRawAfterSuccess','WorkerCount','WorkerIndex', ...
    'OverviewRefreshSeconds','OverviewMaxFiles','ThreadsPerWorker','DetectionMode','TxReferenceChannel'};
values={false,10,20,'auto','',true,false,1,1,60,120,1,'all',[]};
for j=1:numel(names), if ~isfield(config,names{j}),config.(names{j})=values{j};end,end
assert(config.WorkerCount>=1 && config.WorkerCount==fix(config.WorkerCount),'Invalid WorkerCount');
assert(any(strcmp(config.DetectionMode,{'all','transmit'})),'Invalid DetectionMode');
assert(config.WorkerIndex>=1 && config.WorkerIndex<=config.WorkerCount && ...
    config.WorkerIndex==fix(config.WorkerIndex),'Invalid WorkerIndex');
assert(exist(input_folder,'dir')==7,'Input folder does not exist');
ensure_directory(output_folder);
if exist('OCTAVE_VERSION','builtin'), fftw('threads',config.ThreadsPerWorker);
elseif exist('maxNumCompThreads','builtin') || exist('maxNumCompThreads','file')
    maxNumCompThreads(config.ThreadsPerWorker);
end
previous=struct('name',{},'bytes',{},'datenum',{});
stats=struct('processed',0,'skipped_complete',0,'deleted_raw',0,'failed',0,'other_worker',0);
started=tic; last_overview=-Inf;
while true
    files=dir(fullfile(input_folder,'20*.raw'));
    [~,order]=sort({files.name}); files=files(order);
    for j=1:numel(files)
        file=files(j);
        if source_worker(file.name,config.WorkerCount)~=config.WorkerIndex
            stats.other_worker=stats.other_worker+1; continue;
        end
        [~,stem]=fileparts(file.name); raw_path=fullfile(input_folder,file.name);
        try
            day=datestr(datenum(stem(1:14),'yyyymmddHHMMSS')-8/24,'yyyy-mm-dd');
        catch problem
            stats.failed=stats.failed+1; warning('Invalid source date: %s',problem.message); continue;
        end
        day_folder=fullfile(output_folder,day); target=fullfile(day_folder,stem);
        legacy=fullfile(output_folder,stem);
        if exist(fullfile(legacy,'catalogue.mat'),'file')==2 && ~exist(target,'dir'), target=legacy; end
        done=fullfile(target,'catalogue.mat');
        if exist(done,'file')==2
            record=load(done,'result'); r=record.result;
            stored_mode='all';
            if isfield(r,'detection_mode'), stored_mode=r.detection_mode; end
            if ~strcmp(stored_mode,config.DetectionMode)
                warning('Detection mode changed: use a new output folder; raw retained');
                stats.failed=stats.failed+1; continue;
            end
            if r.complete && r.source_bytes==file.bytes && r.source_datenum==file.datenum
                stats.skipped_complete=stats.skipped_complete+1;
                if config.MakeOverview && isfield(r,'overview') && ...
                        (~r.overview_written || exist(fullfile(target,'overview.png'),'file')~=2 || ...
                        exist(fullfile(target,'overview.pdf'),'file')~=2)
                    try, r=repair_overview(target,r);
                    catch problem, warning('Overview repair failed; raw retained: %s',problem.message); end
                end
                if config.DeleteRawAfterSuccess
                    try, stats.deleted_raw=stats.deleted_raw+cleanup_raw(raw_path,target,r,config);
                    catch problem, warning('Raw retained: %s',problem.message); end
                end
            else, warning('Source changed or incomplete result: %s; results preserved',file.name); end
            continue;
        end
        if exist(target,'dir')==7 || exist([target '.partial'],'dir')==7
            warning('Incomplete output preserved for inspection: %s',target); continue;
        end
        old=find(strcmp({previous.name},file.name),1);
        stable=~isempty(old) && previous(old).bytes==file.bytes && previous(old).datenum==file.datenum;
        if config.RunOnce, stable=true; end
        if ~stable || (now-file.datenum)*86400<config.StableSeconds, continue; end
        if ~config.RunOnce && (now-file.datenum)*86400>180
            warning('Processing backlog: %s has waited more than 180 seconds',file.name);
        end
        ensure_directory(day_folder);
        claim=[target '.lock']; [ok,message]=mkdir(claim);
        if ~ok || ~isempty(message), continue; end
        release=onCleanup(@()release_claim(claim));
        options=struct('Format',config.Format,'MakeOverview',config.MakeOverview, ...
            'DetectionMode',config.DetectionMode,'TxReferenceChannel',config.TxReferenceChannel);
        if ~isempty(config.PointingFolder)
            options.PointingFile=fullfile(config.PointingFolder,[stem '_pointing.mat']);
        end
        try
            r=tibet_process_file(raw_path,target,options);
            stats.processed=stats.processed+1;
            if config.DeleteRawAfterSuccess
                stats.deleted_raw=stats.deleted_raw+cleanup_raw(raw_path,target,r,config);
            end
        catch problem
            stats.failed=stats.failed+1;
            warning('Raw retained; processing/cleanup failed: %s: %s',file.name,problem.message);
        end
        clear release;
    end
    if config.MakeOverview && config.WorkerIndex==1 && ...
            (config.RunOnce || toc(started)-last_overview>=config.OverviewRefreshSeconds)
        try, tibet_refresh_overviews(output_folder,config.OverviewMaxFiles);
        catch problem, warning('Daily overview failed: %s',problem.message); end
        last_overview=toc(started);
    end
    if config.RunOnce, break; end
    previous=files; pause(config.PollSeconds);
end
stats.elapsed_seconds=toc(started);
end

function worker=source_worker(name,count)
[~,stem]=fileparts(name); token=regexp(stem,'_([0-9]+)(?:\([0-9]+\))?$','tokens','once');
if ~isempty(token), number=str2double(token{1});
else
    number=0;
    for c=double(stem), number=mod(number*131+c,2147483647); end
end
worker=mod(number,count)+1;
end

function ensure_directory(path)
if exist(path,'dir')~=7
    [ok,message]=mkdir(path);
    assert(ok || exist(path,'dir')==7,message);
end
end

function release_claim(path)
if exist(path,'dir')==7, rmdir(path); end
end

function deleted=cleanup_raw(path,target,result,config)
deleted=0;
assert(result.complete && isfield(result,'cut_storage_verified') && result.cut_storage_verified, ...
    'No verified completed cut set; refusing source deletion');
assert(~config.MakeOverview || result.overview_written,'Overview missing; refusing source deletion');
if config.MakeOverview
    assert(exist(fullfile(target,'overview.png'),'file')==2 && ...
        exist(fullfile(target,'overview.pdf'),'file')==2,'Overview files missing; refusing deletion');
end
assert(numel(result.event_files)==result.event_count,'Incomplete event list');
assert(isfield(result,'event_checksums') && numel(result.event_checksums)==result.event_count, ...
    'No saved event-integrity checksums; refusing deletion');
for j=1:numel(result.event_files)
    item=dir(fullfile(target,result.event_files{j}));
    assert(numel(item)==1 && item.bytes>0,'A verified event file is missing');
    assert(tibet_file_checksum(fullfile(target,result.event_files{j}))==result.event_checksums(j), ...
        'Saved event changed; refusing source deletion');
end
item=dir(path);
assert(numel(item)==1 && item.bytes==result.source_bytes && item.datenum==result.source_datenum, ...
    'Source changed after processing; refusing deletion');
fid=fopen(path,'rb','ieee-le'); assert(fid>=0,'Cannot reopen source');
header=fread(fid,8192,'*uint8'); fclose(fid);
assert(isequal(header,result.header),'Source header changed; refusing deletion');
retention=struct('source_path',path,'source_bytes',result.source_bytes, ...
    'event_files',{result.event_files},'event_count',result.event_count, ...
    'cut_storage_verified',true,'raw_deleted',false,'recorded_local_time',datestr(now,'yyyy-mm-dd HH:MM:SS'));
save_retention(fullfile(target,'retention.mat'),retention);
delete(path); assert(exist(path,'file')~=2,'Source deletion failed');
retention.raw_deleted=true; save_retention(fullfile(target,'retention.mat'),retention);
deleted=1;
end

function save_retention(path,retention)
if exist('OCTAVE_VERSION','builtin'), save('-mat7-binary',path,'retention');
else, save(path,'retention','-v7'); end
end

function result=repair_overview(target,result)
claim=[target '.lock']; [ok,message]=mkdir(claim);
assert(ok && isempty(message),'Another worker owns overview repair');
release=onCleanup(@()release_claim(claim));
tibet_plot_overview({result},fullfile(target,'overview'));
result.overview_written=true;
temporary=fullfile(target,'catalogue_repaired.mat');
if exist('OCTAVE_VERSION','builtin'), save('-mat7-binary',temporary,'result');
else, save(temporary,'result','-v7'); end
[ok,message]=movefile(temporary,fullfile(target,'catalogue.mat'),'f'); assert(ok,message);
end
