function tibet_plot_overview(results,stem)
% Compact range-time matched-filter overview and retained detection tracks.
assert(~isempty(results),'No catalogues to plot');
start=min(cellfun(@(r)double(r.filename_epoch_us)/1e6,results));
fig=figure('Visible','off','Color','w','Position',[50 50 1200 900]);
closer=onCleanup(@()close(fig));
set(fig,'PaperUnits','centimeters','PaperPosition',[0 0 18 14], ...
    'PaperSize',[18 14],'InvertHardcopy','off');
upper=axes('Parent',fig,'Position',[0.11 0.57 0.75 0.36]); hold(upper,'on');
lower=axes('Parent',fig,'Position',[0.11 0.11 0.75 0.32]); hold(lower,'on');
colors=lines(7); number=0; stop=0; low=Inf; high=-Inf; lower_min=Inf; lower_max=-Inf;
for j=1:numel(results)
    r=results{j}; offset=double(r.filename_epoch_us)/1e6-start;
    t=offset+double(r.overview.pulse_index)*r.ipp_us*1e-6;
    db=10*log10(max(double(r.overview.peak_power_ratio),1e-12));
    imagesc(upper,t,r.overview.range_km,db);
    half=median(diff(r.overview.range_km))/2;
    low=min(low,min(r.overview.range_km)-half); high=max(high,max(r.overview.range_km)+half);
    stop=max(stop,offset+r.processed_pulses*r.ipp_us*1e-6);
    for k=1:numel(r.events)
        p=r.events{k}; number=number+1;
        x=offset+[p.pulse_index]*r.ipp_us*1e-6;
        range=[p.delay_samples]*299792458/(2*r.sample_rate_hz*1000);
        lower_min=min(lower_min,min(range)); lower_max=max(lower_max,max(range));
        plot(lower,x,range,'.','Color',colors(mod(number-1,7)+1,:),'MarkerSize',7);
    end
end
set([upper lower],'FontSize',11,'FontName','Times','Box','on','Layer','top');
set(upper,'YDir','normal'); caxis(upper,[0 25]);
ylim(upper,[low high]);
if exist('parula','file') || exist('parula','builtin'), colormap(fig,parula(256));
else, colormap(fig,viridis(256)); end
bar=colorbar('peer',upper); set(bar,'FontName','Times','FontSize',11);
ylabel(bar,'Peak / median power (dB)','FontSize',11,'FontName','Times');
% Colorbar creation must not resize the aligned time axes.
set(upper,'Position',[0.11 0.57 0.75 0.36]);
set(bar,'Position',[0.89 0.57 0.018 0.36]);
set([upper lower],'XLim',[0 max(stop,0.001)]);
ylabel(upper,'Nominal slant range (km)'); ylabel(lower,'Nominal slant range (km)');
xlabel(upper,'Time since start (s)'); xlabel(lower,'Time since start (s)');
title(upper,[datestr(datenum(1970,1,1)+start/86400,'yyyy-mm-dd HH:MM:SS') ' UTC'], ...
    'FontSize',11,'FontWeight','normal');
if number==1, label='1 retained candidate track'; else, label=sprintf('%d retained candidate tracks',number); end
title(lower,label,'FontSize',11,'FontWeight','normal');
grid(lower,'on');
if number==0, ylim(lower,[70 180]);
else, margin=max(0.3,0.05*(lower_max-lower_min)); ylim(lower,[lower_min-margin lower_max+margin]); end
temporary=[stem '.partial'];
drawnow();
print(fig,[temporary '.pdf'],'-dpdf','-painters');
print(fig,[temporary '.png'],'-dpng','-r160','-painters');
[ok,message]=movefile([temporary '.png'],[stem '.png'],'f'); assert(ok,message);
[ok,message]=movefile([temporary '.pdf'],[stem '.pdf'],'f'); assert(ok,message);
end
