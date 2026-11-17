function tibet_verify_saved_event(path,event)
% Reopen the actual saved voltage and essential metadata before permitting cleanup.
[~,~,extension]=fileparts(path);
is_mat=strcmpi(extension,'.mat');
if is_mat, stored=load(path); end
fields={'zrx_echoes_iq_i16','zrx_echoes_overflow_index','zrx_echoes_overflow_words_i16', ...
    'ztx_history_iq_i16','ztx_history_overflow_index','ztx_history_overflow_words_i16', ...
    'pulse_index','record_counter','ipp_start_time_us','raw_header','source_channel_indices', ...
    'sample_rate_hz','frequency_hz','range_start_gate','tx_start_gate','tx_history_valid', ...
    'tx_history_decoded_baud_sign','tx_alpha_per_pulse','tx_beta_per_pulse','tx_beam_id'};
for j=1:numel(fields)
    key=fields{j}; expected=event.(key);
    if is_mat, actual=stored.(key);
    elseif isempty(expected)
        assert(h5read(path,['/' key '_empty'])==1,'Missing empty-array marker'); actual=expected;
    else, actual=h5read(path,['/' key]); end
    assert(isequal(actual(:),expected(:)),'Saved event verification failed: %s',key);
end
for key={'pulse_index','delay_samples','doppler_hz','matched_amplitude_re_im','noise_power','tx_time_us'}
    field=key{1}; expected=event.detections.(field);
    if is_mat, actual=stored.detections.(field); else, actual=h5read(path,['/detections/' field]); end
    assert(isequal(actual(:),expected(:)),'Saved detection verification failed: %s',field);
end
end
