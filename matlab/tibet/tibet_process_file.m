function result = tibet_process_file(raw_path, output_dir, options)
% Standalone MATLAB/Octave range-Doppler meteor detector. No toolboxes required.
if nargin < 3, options = struct(); end
options = defaults(options);
if exist(output_dir, 'dir'), error('Output already exists: %s', output_dir); end
fid = fopen(raw_path, 'rb', 'ieee-le');
if fid < 0, error('Cannot open %s', raw_path); end
closer = onCleanup(@() fclose(fid));
header = fread(fid, 8192, '*uint8');
h = parse_header(header, raw_path);
info = dir(raw_path);
assert(info.bytes >= 8192+h.pulses*h.record_bytes, 'Incomplete declared records');
assert(h.beam_total <= 1 || ~isempty(options.PointingFile), ...
    'Multi-beam mode requires a per-pulse PointingFile; direction must not be guessed');
pointing = load_pointing(h, options.PointingFile);
gates = 241:608;
windows = bsxfun(@plus, (0:31)', gates+1);
empty = struct('pulse_index', {}, 'delay_samples', {}, 'doppler_bin', {}, ...
    'doppler_hz', {}, 'statistic_db', {}, 'amplitude', {}, 'noise', {}, 'beam', {});
points = empty; tracks = {}; last = []; candidate_count = 0; veto_count = 0;
overview=struct('pulse_index',uint64(0:16:min(h.pulses,options.MaxPulses)-1), ...
    'range_km',(gates-7)*299792458/(2*h.fs*1000), ...
    'peak_power_ratio',zeros(numel(gates),ceil(min(h.pulses,options.MaxPulses)/16),'single'));
tic;
for batch = 0:options.BatchPulses:min(h.pulses,options.MaxPulses)-1
    count = min([options.BatchPulses, h.pulses-batch, options.MaxPulses-batch]);
    raw = read_records(fid, h, batch, count);
    for bi = 1:count
        pulse = batch+bi-1;
        re = reshape(raw(1,:,:,bi), 15, 640).';
        im = reshape(raw(2,:,:,bi), 15, 640).';
        z = complex(single(re),single(im));
        noise_power = sort(sq(z(251:640,:)),1);
        noise = max(noise_power(196,:)/single(log(2)), single(1));
        tx = z(8:39,3); energy = sum(sq(tx));
        assert(energy > 3*32*noise(3), 'Weak transmit reference at pulse %d', pulse);
        tx = tx/sqrt(energy);
        rx = bsxfun(@rdivide,z(:,1:5),sqrt(noise(1:5)));
        segments = reshape(rx(windows(:),:),32,numel(gates),5);
        spectrum = fftshift(fft(bsxfun(@times,segments,conj(tx)),128,1),1);
        power = sum(sq(spectrum),3); reference = median(power(:));
        if mod(pulse,16)==0, overview.peak_power_ratio(:,pulse/16+1)=max(power,[],1)'/reference; end
        indices = find(power/reference-1 >= single(10^(options.ThresholdDB/10)));
        [bin,row] = ind2sub(size(power),indices);
        rr = max(1,min(numel(gates),bsxfun(@plus,row(:),[-1 -1 -1 0 0 0 1 1 1])));
        kk = mod(bsxfun(@plus,bin(:)-1,[-1 0 1 -1 0 1 -1 0 1]),128)+1;
        other = kk+128*(rr-1);
        local = all(bsxfun(@lt,power(other),power(indices)) | ...
            (bsxfun(@eq,power(other),power(indices)) & bsxfun(@ge,other,indices)),2);
        indices = indices(local);
        if ~isempty(indices)
            [~,order] = sortrows([-double(power(indices)),double(indices)],[1 2]);
            indices = indices(order);
        end
        selected = zeros(0,2);
        for j = 1:numel(indices)
            [bin,row] = ind2sub(size(power),indices(j));
            delta_bin = abs(selected(:,1)-bin);
            if any(abs(selected(:,2)-row)<=5 & min(delta_bin,128-delta_bin)<=8), continue; end
            selected(end+1,:) = [bin,row]; %#ok<AGROW>
            if size(selected,1)==3, break; end
        end
        accepted = empty;
        for j = 1:size(selected,1)
            bin = selected(j,1); row = selected(j,2);
            logp = double(log(max(power(mod((bin-2):bin,128)+1,row),single(1e-20))));
            curvature = logp(1)-2*logp(2)+logp(3); delta = 0;
            if curvature < 0, delta = max(-0.5,min(0.5,(logp(1)-logp(3))/(2*curvature))); end
            frequency = (bin-1-64+delta)*h.fs/128;
            phase = -single(2*pi)*single(frequency)*single((0:31)')/single(h.fs);
            phasor = complex(cos(phase),sin(phase));
            amplitude = sum(bsxfun(@times,z(gates(row)+(1:32),:),conj(tx).*phasor),1);
            [~,beam] = max(sq(amplitude(1:5))./noise(1:5));
            p = struct('pulse_index',pulse,'delay_samples',gates(row)-7, ...
                'doppler_bin',bin-1,'doppler_hz',frequency, ...
                'statistic_db',10*log10(power(bin,row)/reference-1), ...
                'amplitude',amplitude,'noise',noise,'beam',beam-1);
            veto = false;
            for k = 1:numel(accepted)
                if code_alias(p,accepted(k),tx,h.fs), veto = true; break; end
            end
            if veto, veto_count = veto_count+1; continue; end
            accepted(end+1)=p; candidate_count=candidate_count+1; %#ok<AGROW>
            points(end+1)=p; point_index=numel(points); %#ok<AGROW>
            best = 0; best_score = Inf;
            if ~isempty(last)
                previous = points(last);
                dt = (pulse-[previous.pulse_index])*h.ipp_us*1e-6;
                error = abs(p.delay_samples*299792458/(2*h.fs*1000) - ...
                    ([previous.delay_samples]*299792458/(2*h.fs*1000)- ...
                    [previous.doppler_hz]*h.wavelength/2000.*dt));
                dv = abs(p.doppler_hz-[previous.doppler_hz])*h.wavelength/2000;
                eligible = find(dt>0 & dt<=0.05 & error<=0.75 & dv<=25);
                for k = eligible
                    prior = previous(k);
                    if any(pointing.alpha(pulse+1)~=pointing.alpha(prior.pulse_index+1)) || ...
                            any(pointing.beta(pulse+1)~=pointing.beta(prior.pulse_index+1)) || ...
                            pointing.beam_id(pulse+1)~=pointing.beam_id(prior.pulse_index+1), continue; end
                    a=p.amplitude(9:15)./sqrt(p.noise(9:15));
                    b=prior.amplitude(9:15)./sqrt(prior.noise(9:15));
                    coherence=abs(sum(a.*conj(b)))/max(sqrt(sum(sq(a))*sum(sq(b))),single(1e-20));
                    score=error(k)+dv(k)/50+1-double(coherence);
                    if coherence>=0.65 && score<best_score, best=k; best_score=score; end
                end
            end
            if best==0, tracks{end+1}=point_index; last(end+1)=point_index; %#ok<AGROW>
            else, tracks{best}(end+1)=point_index; last(best)=point_index; end
        end
    end
    if options.Verbose && mod(batch,512)==0
        fprintf('Pulse %d/%d; %d candidates; %.1f s\n',batch+count,h.pulses,candidate_count,toc);
    end
end
events = {};
for j=1:numel(tracks)
    p=points(tracks{j});
    if numel(p)<8 || (p(end).pulse_index-p(1).pulse_index)*h.ipp_us<30000, continue; end
    t=([p.pulse_index]-p(1).pulse_index)*h.ipp_us*1e-6;
    range=[p.delay_samples]*299792458/(2*h.fs*1000);
    rate=sum((t-mean(t)).*(range-mean(range)))/sum((t-mean(t)).^2);
    doppler_velocity=mean([p.doppler_hz])*h.wavelength/2000;
    if abs(rate+doppler_velocity)<10, events{end+1}=p; end %#ok<AGROW>
end
partial=[output_dir '.partial'];
assert(~exist(partial,'dir'),'Incomplete output exists; inspect it before retrying');
mkdir(partial);
event_files=cell(1,numel(events)); event_bytes=0; event_checksums=zeros(1,numel(events),'uint32');
for j=1:numel(events)
    event=make_cut(fid,h,pointing,events{j},j-1);
    path=tibet_save_event(fullfile(partial,sprintf('event_%05d',j-1)),event,options.Format);
    tibet_verify_saved_event(path,event);
    event_checksums(j)=tibet_file_checksum(path);
    [~,base,extension]=fileparts(path); event_files{j}=[base extension];
    item=dir(path); event_bytes=event_bytes+item.bytes;
end
result=struct('schema','tibet_matlab_detector_v2','source_path',raw_path, ...
    'source_bytes',info.bytes,'source_datenum',info.datenum,'header',header, ...
    'candidate_count',candidate_count,'code_alias_veto_count',veto_count, ...
    'event_count',numel(events),'events',{events},'threshold_db',options.ThresholdDB, ...
    'processed_pulses',min(h.pulses,options.MaxPulses),'processing_seconds',toc, ...
    'sample_rate_hz',h.fs,'ipp_us',h.ipp_us,'filename_epoch_us',h.epoch_us, ...
    'overview',overview,'overview_written',false,'event_files',{event_files}, ...
    'event_bytes',event_bytes,'event_checksums',event_checksums,'cut_storage_verified',true, ...
    'complete',min(h.pulses,options.MaxPulses)==h.pulses);
if options.MakeOverview
    try
        tibet_plot_overview({result},fullfile(partial,'overview'));
        result.overview_written=true;
    catch problem
        result.overview_plot_error=problem.message;
        warning('Cuts saved; overview failed: %s',problem.message);
    end
end
result.processing_seconds=toc;
save_mat(fullfile(partial,'catalogue.mat'),'result',result);
[ok,message]=movefile(partial,output_dir); assert(ok,message);
fprintf('%s: %d events, %d points, %.1f s\n',raw_path,numel(events), ...
    sum(cellfun(@numel,events)),result.processing_seconds);
end

function o=defaults(o)
names={'BatchPulses','MaxPulses','ThresholdDB','Verbose','Format','PointingFile','MakeOverview'};
values={16,Inf,12,true,'auto','',true};
for j=1:numel(names), if ~isfield(o,names{j}), o.(names{j})=values{j}; end, end
end

function h=parse_header(b,path)
assert(numel(b)==8192 && word(b,0,'uint32')==hex2dec('aaaa5555') && ...
    word(b,8188,'uint32')==hex2dec('aa5555aa'),'Bad raw header');
h=struct('header',b,'path',path,'pulses',word(b,80,'uint32'), ...
    'gates',word(b,240,'uint32'),'channels',word(b,426,'uint16'), ...
    'ipp_us',word(b,76,'uint32'),'frequency',word(b,84,'uint32')*120e6/2^32, ...
    'beam_total',word(b,42,'uint16'),'frame',word(b,44,'uint16'), ...
    'alpha',word(b,256,'uint16'),'beta',word(b,258,'uint16'));
h.fs=8e6/2^double(b(95)); h.wavelength=299792458/h.frequency;
assert(h.fs==500000 && h.gates==640 && h.channels==15 && h.ipp_us==1280 && ...
    h.pulses>0 && abs(h.frequency-50e6)<1000 && word(b,730,'uint16')==4 && ...
    word(b,88,'uint16')==3 && word(b,90,'uint32')==80 && word(b,236,'uint32')==0, ...
    'Unvalidated radar mode: do not silently apply this detector');
assert(isequal(double(typecast(b(549:578),'uint16'))', ...
    [15 15 15 15 15 23 24 25 26 27 28 29 30 31 32]),'Changed channel mapping');
h.record_bytes=640*15*8+16;
[~,name]=fileparts(path); assert(numel(name)>=14,'Missing filename date');
date=datenum(name(1:14),'yyyymmddHHMMSS');
assert(strcmp(datestr(date,'yyyymmddHHMMSS'),name(1:14)),'Invalid filename calendar');
h.epoch_us=uint64(round((date-datenum(1970,1,1))*86400-8*3600))*uint64(1000000);
end

function value=word(b,offset,type)
if strcmp(type,'uint16'), n=2; else, n=4; end
value=double(typecast(b(offset+(1:n)),type));
[~,~,endian]=computer; if endian=='B', value=double(swapbytes(cast(value,type))); end
end

function iq=read_records(fid,h,first,count)
assert(fseek(fid,8192+first*h.record_bytes,'bof')==0,'Seek failed');
words=fread(fid,[h.record_bytes/4,count],'*int32');
assert(size(words,2)==count,'Short pulse record');
magic=typecast(uint32(hex2dec('aabbccdd')),'int32');
assert(all(words(1,:)==magic) && all(words(end-1,:)==magic), 'Bad pulse markers');
assert(all(double(words(2,:))==first+(1:count)) && ...
    all(double(words(end,:))==first+(1:count)),'Pulse counter mismatch');
iq=reshape(words(3:end-2,:),2,15,640,count);
end

function p=sq(z)
p=real(z).^2+imag(z).^2;
end

function yes=code_alias(weak,strong,tx,fs)
shift=weak.delay_samples-strong.delay_samples;
if abs(shift)>=32, yes=false; return; end
k=(0:31)'; j=k+shift; valid=j>=0 & j<32; k=k(valid); j=j(valid);
phase=single(2*pi*(strong.doppler_hz*j-weak.doppler_hz*k)/fs);
response=sum(tx(j+1).*conj(tx(k+1)).*complex(cos(phase),sin(phase)));
observed=sum(sq(weak.amplitude)./weak.noise);
residual=sum(sq(weak.amplitude-strong.amplitude*response)./weak.noise);
yes=observed>0 && residual/observed<0.35;
end

function p=load_pointing(h,file)
p=struct('alpha',repmat(h.alpha,h.pulses,1),'beta',repmat(h.beta,h.pulses,1), ...
    'beam_id',zeros(h.pulses,1,'int32'),'zenith_deg',NaN(h.pulses,1), ...
    'azimuth_deg',NaN(h.pulses,1));
if h.alpha==0 && h.beta==0, p.zenith_deg(:)=0; end
if ~isempty(file)
    a=load(file); assert(isfield(a,'pointing'),'Pointing file must contain pointing struct');
    fields=fieldnames(p);
    for j=1:numel(fields)
        value=a.pointing.(fields{j}); assert(numel(value)==h.pulses,'Pointing schedule length');
        p.(fields{j})=value(:);
    end
    assert(all(isfinite(p.zenith_deg) & p.zenith_deg>=0 & p.zenith_deg<=90),'Invalid zenith angles');
    assert(all(p.zenith_deg==0 | (isfinite(p.azimuth_deg) & p.azimuth_deg>=0 & p.azimuth_deg<360)),'Invalid azimuth');
end
end

function event=make_cut(fid,h,pointing,p,id)
pad=ceil(80000/h.ipp_us); first=max(0,p(1).pulse_index-pad);
last=min(h.pulses-1,p(end).pulse_index+pad);
lo=max(0,min([p.delay_samples])+7-64); hi=min(640,max([p.delay_samples])+7+32+64);
n=last-first+1; width=hi-lo;
raw=read_records(fid,h,first,n);
iq=permute(raw(:,6:15,lo+1:hi,:),[1 3 2 4]);
history=zeros(2,32,3,n,'int32'); valid=zeros(3,n,'uint8');
origin=max(0,first-2); prior=read_records(fid,h,origin,last-origin+1);
tx_alpha=zeros(3,n,'uint16'); tx_beta=tx_alpha; beam=-ones(3,n,'int32');
zenith=NaN(3,n); azimuth=zenith; code=zeros(8,3,n,'int8');
for i=1:n
    pulse=first+i-1;
    for order=0:2
        from=pulse-order; if from<0, continue; end
        history(:,:,order+1,i)=reshape(prior(:,3,8:39,from-origin+1),2,32);
        valid(order+1,i)=1; tx_alpha(order+1,i)=pointing.alpha(from+1);
        tx_beta(order+1,i)=pointing.beta(from+1); beam(order+1,i)=pointing.beam_id(from+1);
        zenith(order+1,i)=pointing.zenith_deg(from+1); azimuth(order+1,i)=pointing.azimuth_deg(from+1);
        z=complex(double(history(1,:,order+1,i)),double(history(2,:,order+1,i)));
        baud=sum(reshape(z,4,8),1); code(:,order+1,i)=int8(2*(real(baud*conj(baud(1)))>=0)-1);
    end
end
[compact,index,words]=tibet_pack_iq(iq);
[tx_compact,tx_index,tx_words]=tibet_pack_iq(history);
event=struct('schema','tibet_matlab_cut_v1','event_id',sprintf('event_%05d',id), ...
    'source_path',h.path,'raw_header',h.header,'source_channel_indices',uint16(5:14), ...
    'channel_names','S0 S1 S2 L0 L1 L2 L3 L4 L5 L6','tx_reference_source_channel',uint16(2), ...
    'sample_rate_hz',uint32(h.fs),'frequency_hz',h.frequency,'ipp_us',uint32(h.ipp_us), ...
    'range_start_gate',uint32(lo),'tx_start_gate',uint32(7),'txlen',uint32(32), ...
    'samples_per_baud',uint32(4),'baud_count',uint32(8),'code_duration_us',64, ...
    'pulse_index',uint64(first:last),'record_counter',uint32(first+1:last+1), ...
    'ipp_start_time_us',h.epoch_us+uint64(first:last)*uint64(h.ipp_us), ...
    'tx_alpha_per_pulse',pointing.alpha(first+1:last+1),'tx_beta_per_pulse',pointing.beta(first+1:last+1), ...
    'tx_beam_id',pointing.beam_id(first+1:last+1),'tx_zenith_angle_deg',pointing.zenith_deg(first+1:last+1), ...
    'tx_azimuth_deg',pointing.azimuth_deg(first+1:last+1), ...
    'tx_history_valid',valid,'tx_history_alpha_code',tx_alpha,'tx_history_beta_code',tx_beta, ...
    'tx_history_beam_id',beam,'tx_history_zenith_angle_deg',zenith,'tx_history_azimuth_deg',azimuth, ...
    'tx_history_decoded_baud_sign',code,'zrx_echoes_iq_i16',compact, ...
    'zrx_echoes_overflow_index',index,'zrx_echoes_overflow_words_i16',words, ...
    'ztx_history_iq_i16',tx_compact,'ztx_history_overflow_index',tx_index, ...
    'ztx_history_overflow_words_i16',tx_words,'detections',point_arrays(p), ...
    'time_status','Filename assumed Beijing UTC+8; binary calendar invalid; absolute time precision unconfirmed', ...
    'phase_status','No receiver phase correction applied; original integer voltages preserved', ...
    'range_alias_status','Previous two codes sampled; even-IPP aliases unresolved', ...
    'iq_layout','MATLAB [I/Q,gate,channel,pulse]; overflow indices zero-based linear, low/high int16 words', ...
    'latitude_deg',30.1,'longitude_deg',90.5,'altitude_m',4287);
event.detections.tx_time_us=h.epoch_us+uint64([p.pulse_index])*uint64(h.ipp_us);
alias_power=NaN(3,numel(p)); replay=zeros(2,10,numel(p),'single');
for j=1:numel(p)
    row=p(j).pulse_index-first+1; gate=p(j).delay_samples+7-lo;
    echo=complex(single(reshape(iq(1,gate+(1:32),:,row),32,10)), ...
        single(reshape(iq(2,gate+(1:32),:,row),32,10)));
    phase=-single(2*pi)*single(p(j).doppler_hz)*single((0:31)')/single(h.fs);
    phasor=complex(cos(phase),sin(phase));
    for order=1:3
        if ~valid(order,row), continue; end
        tx=complex(single(history(1,:,order,row)),single(history(2,:,order,row))).';
        tx=tx/sqrt(sum(sq(tx)));
        amplitude=sum(bsxfun(@times,echo,conj(tx).*phasor),1);
        alias_power(order,j)=sum(double(sq(amplitude))./double(p(j).noise(6:15)));
        if order==1
            replay(:,:,j)=[real(amplitude);imag(amplitude)];
            assert(max(abs(amplitude-p(j).amplitude(6:15))./max(abs(p(j).amplitude(6:15)),1))<1e-5, ...
                'Cut does not reproduce detector amplitudes');
        end
    end
end
event.detections.cut_replayed_amplitude_re_im=replay;
event.detections.alias_matched_power_at_observed_delay_doppler=alias_power;
event.detections.alias_order=uint32(0:2);
event.detections.range_alias_interval_km=299792458*h.ipp_us*1e-6/2000;
geometry_file=fullfile(fileparts(mfilename('fullpath')),'array_geometry.mat');
if exist(geometry_file,'file')==2, event.geometry=load(geometry_file); end
end

function d=point_arrays(p)
d=struct('pulse_index',uint64([p.pulse_index]),'delay_samples',uint32([p.delay_samples]), ...
    'doppler_bin',uint32([p.doppler_bin]),'doppler_hz',[p.doppler_hz], ...
    'detector_excess_over_median_db',single([p.statistic_db]), ...
    'matched_amplitude_re_im',permute(cat(3,real(vertcat(p.amplitude)),imag(vertcat(p.amplitude))),[3 2 1]), ...
    'noise_power',vertcat(p.noise).','strongest_full_array_rx',uint16([p.beam]));
d.full_array_matched_amplitude_re_im=d.matched_amplitude_re_im(:,1:5,:);
d.full_array_noise_power=d.noise_power(1:5,:);
d.matched_amplitude_re_im=d.matched_amplitude_re_im(:,6:15,:);
d.noise_power=d.noise_power(6:15,:);
end

function save_mat(path,name,value)
s=struct(); s.(name)=value;
if exist('OCTAVE_VERSION','builtin'), save('-mat7-binary',path,'-struct','s');
else, save(path,'-struct','s','-v7'); end
end
