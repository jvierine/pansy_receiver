function checksum=tibet_file_checksum(path)
% Streaming Adler-32 detects accidental changes after verified event publication.
fid=fopen(path,'rb'); assert(fid>=0,'Cannot checksum %s',path);
closer=onCleanup(@()fclose(fid)); a=1; b=0;
while true
    bytes=fread(fid,65536,'*uint8'); n=numel(bytes); if n==0, break; end
    values=double(bytes(:));
    b=mod(b+n*a+sum((n:-1:1)'.*values),65521);
    a=mod(a+sum(values),65521);
end
checksum=uint32(b*65536+a);
end
