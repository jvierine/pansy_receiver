function path=tibet_save_event(stem,event,format)
% MATLAB HDF5 uses shuffle then deflate. Octave fallback is compressed MAT v7.
if strcmp(format,'auto')
    if exist('h5create','file')==2, format='hdf5'; else, format='mat'; end
end
if strcmp(format,'hdf5')
    assert(exist('h5create','file')==2,'HDF5 writing requires h5create');
    path=[stem '.h5']; write_group(path,'',event);
elseif strcmp(format,'mat')
    path=[stem '.mat'];
    if exist('OCTAVE_VERSION','builtin'), save('-mat7-binary',path,'-struct','event');
    else, save(path,'-struct','event','-v7'); end
else, error('Format must be auto, hdf5 or mat'); end
end

function write_group(path,prefix,s)
fields=fieldnames(s);
for i=1:numel(fields)
    name=fields{i}; value=s.(name); key=[prefix '/' name];
    if isstruct(value), write_group(path,key,value); continue; end
    if ischar(value), value=uint8(value); end
    if isempty(value)
        % Empty overflow tables are identified by this scalar count in HDF5.
        h5create(path,[key '_empty'],1,'Datatype','uint8'); h5write(path,[key '_empty'],uint8(1));
        continue;
    end
    if islogical(value), value=uint8(value); end
    shape=size(value); chunk=shape; chunk(end)=min(chunk(end),16);
    h5create(path,key,shape,'Datatype',class(value),'ChunkSize',chunk, ...
        'Shuffle',true,'Deflate',4);
    h5write(path,key,value);
end
end
