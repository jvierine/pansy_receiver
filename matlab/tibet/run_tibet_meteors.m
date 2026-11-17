% Edit these folders. No external executable, Python, or MATLAB toolbox needed.
InputFolder = 'D:\TibetMST\raw';
OutputFolder = 'D:\TibetMST\meteors';
RunOnce = false;
PollSeconds = 10;

config = struct('RunOnce',RunOnce,'PollSeconds',PollSeconds,'StableSeconds',20, ...
    'Format','auto','PointingFolder','');
tibet_watch(InputFolder,OutputFolder,config);
