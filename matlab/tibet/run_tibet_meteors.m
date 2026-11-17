% Edit these folders. No external executable, Python, or MATLAB toolbox needed.
InputFolder = 'D:\TibetMST\raw';
OutputFolder = 'D:\TibetMST\meteors';
RunOnce = false;
PollSeconds = 10;
if ~exist('WorkerCount','var'), WorkerCount = 1; end
if ~exist('WorkerIndex','var'), WorkerIndex = 1; end
DeleteRawAfterSuccess = false; % Enable only after a successful site pilot.

config = struct('RunOnce',RunOnce,'PollSeconds',PollSeconds,'StableSeconds',20, ...
    'Format','auto','PointingFolder','','MakeOverview',true, ...
    'WorkerCount',WorkerCount,'WorkerIndex',WorkerIndex, ...
    'DeleteRawAfterSuccess',DeleteRawAfterSuccess);
tibet_watch(InputFolder,OutputFolder,config);
