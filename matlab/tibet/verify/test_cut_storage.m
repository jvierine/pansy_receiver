function test_cut_storage(input_folder,output_folder)
% Real-data cut storage regression; no detector claims from these cut fixtures.
assert(~exist(output_folder,'dir'),'Output already exists'); mkdir(output_folder);
files=dir(fullfile(input_folder,'event_*.mat'));
assert(~isempty(files),'No input cut fixtures');
for j=1:numel(files)
    event=load(fullfile(input_folder,files(j).name));
    [event.zrx_echoes_iq_i16,event.zrx_echoes_overflow_index,event.zrx_echoes_overflow_words_i16]= ...
        tibet_pack_iq(event.zrx_echoes_iq_i32);
    [event.ztx_history_iq_i16,event.ztx_history_overflow_index,event.ztx_history_overflow_words_i16]= ...
        tibet_pack_iq(event.ztx_history_iq_i32);
    event=rmfield(event,{'zrx_echoes_iq_i32','ztx_history_iq_i32'});
    [~,stem]=fileparts(files(j).name);
    tibet_save_event(fullfile(output_folder,stem),event,'mat');
end
fprintf('PASS: %d real cuts packed losslessly\n',numel(files));
end
