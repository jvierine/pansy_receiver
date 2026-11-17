function tibet_refresh_overviews(root,max_files)
% Only worker 1 refreshes rolling day plots; per-source products are immutable.
if nargin<2, max_files=120; end
days=dir(root); days=days([days.isdir]);
keep=cellfun(@(s)~isempty(regexp(s,'^[0-9]{4}-[0-9]{2}-[0-9]{2}$','once')),{days.name});
days=days(keep);
[~,order]=sort({days.name}); days=days(order);
for i=max(1,numel(days)-1):numel(days)
    day=fullfile(root,days(i).name); entries=dir(day); entries=entries([entries.isdir]);
    paths={};
    for j=1:numel(entries)
        name=entries(j).name;
        if isempty(name) || name(1)=='.' || ~isempty(regexp(name,'\.(partial|lock)$','once')), continue; end
        path=fullfile(day,name,'catalogue.mat');
        if exist(path,'file')==2, paths{end+1}=path; end %#ok<AGROW>
    end
    paths=sort(paths); paths=paths(max(1,numel(paths)-max_files+1):end);
    results={};
    for j=1:numel(paths)
        record=load(paths{j},'result'); r=record.result;
        if r.complete && isfield(r,'overview'), results{end+1}=r; end %#ok<AGROW>
    end
    if ~isempty(results), tibet_plot_overview(results,fullfile(day,'overview')); end
end
end
