function test_retention(fixture,root)
% Deletes only freshly copied synthetic fixtures, never the original source.
assert(~exist(root,'dir'),'Use a new test directory'); mkdir(root);
cfg=struct('RunOnce',true,'StableSeconds',0,'Format','mat','MakeOverview',false, ...
    'DeleteRawAfterSuccess',true);
input=fullfile(root,'success_input'); mkdir(input); output=fullfile(root,'success_output');
raw=fullfile(input,'20261114060000_fixture.raw'); copyfile(fixture,raw);
s=tibet_watch(input,output,cfg);
assert(s.processed==1 && s.deleted_raw==1 && s.failed==0 && ~exist(raw,'file'));
saved=load(fullfile(output,'2026-11-13','20261114060000_fixture','retention.mat'));
assert(saved.retention.raw_deleted);

input=fullfile(root,'changed_input'); mkdir(input); output=fullfile(root,'changed_output');
raw=fullfile(input,'20261114080000_fixture.raw'); copyfile(fixture,raw);
cfg.DeleteRawAfterSuccess=false; s=tibet_watch(input,output,cfg); assert(s.processed==1);
event=fullfile(output,'2026-11-14','20261114080000_fixture','event_00000.mat');
fid=fopen(event,'ab'); fwrite(fid,uint8(0),'uint8'); fclose(fid);
cfg.DeleteRawAfterSuccess=true; s=tibet_watch(input,output,cfg);
assert(s.deleted_raw==0 && exist(raw,'file')==2,'Corrupted-output safeguard failed');

input=fullfile(root,'empty_input'); mkdir(input); output=fullfile(root,'empty_output');
raw=fullfile(input,'20261114060000_empty.raw'); copyfile(fixture,raw);
fid=fopen(raw,'r+b','ieee-le'); fseek(fid,80,'bof'); fwrite(fid,uint32(16),'uint32'); fclose(fid);
s=tibet_watch(input,output,cfg); assert(s.processed==1 && s.deleted_raw==1 && s.failed==0);
saved=load(fullfile(output,'2026-11-13','20261114060000_empty','catalogue.mat'));
assert(saved.result.event_count==0 && saved.result.complete);

input=fullfile(root,'repair_input'); mkdir(input); output=fullfile(root,'repair_output');
raw=fullfile(input,'20261114060000_repair.raw'); copyfile(fixture,raw);
cfg.DeleteRawAfterSuccess=false; s=tibet_watch(input,output,cfg); assert(s.processed==1);
cfg.MakeOverview=true; cfg.DeleteRawAfterSuccess=true;
s=tibet_watch(input,output,cfg); assert(s.processed==0 && s.deleted_raw==1 && ~exist(raw,'file'));
cfg.MakeOverview=false;

input=fullfile(root,'damaged_input'); mkdir(input); output=fullfile(root,'damaged_output');
raw=fullfile(input,'20261114060000_damaged.raw'); copyfile(fixture,raw);
fid=fopen(raw,'r+b','ieee-le'); fseek(fid,8192,'bof'); fwrite(fid,uint32(0),'uint32'); fclose(fid);
s=tibet_watch(input,output,cfg); assert(s.deleted_raw==0 && s.failed==1 && exist(raw,'file')==2);
for offset=[80 84]
    input=fullfile(root,sprintf('bad_header_%d_input',offset)); mkdir(input);
    output=fullfile(root,sprintf('bad_header_%d_output',offset));
    raw=fullfile(input,'20261114060000_badheader.raw'); copyfile(fixture,raw);
    fid=fopen(raw,'r+b','ieee-le'); fseek(fid,offset,'bof'); fwrite(fid,uint32(0),'uint32'); fclose(fid);
    s=tibet_watch(input,output,cfg);
    assert(s.deleted_raw==0 && s.failed==1 && exist(raw,'file')==2);
end
fprintf('PASS retention: verified deletion, changed-cut refusal, zero-event cleanup, overview repair, damaged input/header retention, UTC dates\n');
end
