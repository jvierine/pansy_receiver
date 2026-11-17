function test_beam_scan(fixture,root)
% Synthetic frame headers exercise the documented scan protocol, not sky data.
assert(~exist(root,'dir'),'Use a new test folder'); mkdir(root);
input=fullfile(root,'raw'); mkdir(input);
fid=fopen(fixture,'rb','ieee-le'); header=fread(fid,8192,'*uint8'); fclose(fid);
alpha=typecast(header(429:438),'uint16'); beta=typecast(header(469:478),'uint16');
for beam=1:5
    path=fullfile(input,sprintf('202611140600%02d_scan_%d.raw',beam,beam));
    copyfile(fixture,path); fid=fopen(path,'r+b','ieee-le');
    fseek(fid,42,'bof'); fwrite(fid,uint16([5 beam-1]),'uint16');
    fseek(fid,256,'bof'); fwrite(fid,uint16([alpha(beam) beta(beam)]),'uint16');
    fclose(fid);
end
output=fullfile(root,'events');
s=tibet_watch(input,output,struct('RunOnce',true,'StableSeconds',0,'MakeOverview',false));
assert(s.processed==5 && s.failed==0);
for beam=1:5
    target=fullfile(output,'2026-11-13',sprintf('202611140600%02d_scan_%d',beam,beam));
    a=load(fullfile(target,'catalogue.mat')); assert(a.result.event_count==1);
    a=load(fullfile(target,'event_00000.mat'));
    assert(all(a.tx_beam_id==beam-1) && all(a.tx_history_beam_id(a.tx_history_valid>0)==beam-1));
    assert(all(a.tx_alpha_per_pulse==alpha(beam)) && all(a.tx_beta_per_pulse==beta(beam)));
    assert(a.tx_reference_source_channel==2 && a.beam_cycle_count==5);
    assert(size(a.zrx_echoes_iq_i16,3)==10);
    path=fullfile(input,sprintf('202611140600%02d_scan_%d.raw',beam,beam));
    a=tibet_process_file(path,fullfile(root,sprintf('matched_%d',beam)), ...
        struct('DetectionMode','transmit','MakeOverview',false,'MaxPulses',16,'Verbose',false));
    assert(all(a.detection_channel==beam));
end
% Real headers show TX/DBF values differing by one steering-code count.
path=fullfile(input,'20261114060004_scan_4.raw'); fid=fopen(path,'r+b','ieee-le');
fseek(fid,256,'bof'); fwrite(fid,uint16(alpha(4)-1),'uint16'); fclose(fid);
a=tibet_process_file(path,fullfile(root,'quantized_match'), ...
    struct('DetectionMode','transmit','MakeOverview',false,'MaxPulses',16,'Verbose',false));
assert(all(a.detection_channel==4));
% An explicit within-frame schedule is still authoritative.
pointing=struct('alpha',repmat(double(alpha(3)),512,1),'beta',repmat(double(beta(3)),512,1), ...
    'beam_id',zeros(512,1),'zenith_deg',zeros(512,1),'azimuth_deg',zeros(512,1));
pointing.alpha(145:end)=alpha(4); pointing.beta(145:end)=beta(4);
pointing.beam_id(145:end)=3; pointing.zenith_deg(145:end)=10; pointing.azimuth_deg(145:end)=90;
schedule=fullfile(root,'pointing.mat');
if exist('OCTAVE_VERSION','builtin'), save('-mat7-binary',schedule,'pointing');
else, save(schedule,'pointing','-v7'); end
a=tibet_process_file(fixture,fullfile(root,'within_frame'), ...
    struct('PointingFile',schedule,'MakeOverview',false,'Verbose',false));
for k=1:numel(a.events)
    pulse=[a.events{k}.pulse_index]+1;
    assert(numel(unique(pointing.beam_id(pulse)))==1,'Track crossed beam change');
end
fprintf('PASS five frame pointings, beam/history IDs, ten-channel cuts, steering tolerance, explicit schedule\n');
end
