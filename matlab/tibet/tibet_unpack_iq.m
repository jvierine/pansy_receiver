function iq=tibet_unpack_iq(compact,index,words)
iq=int32(compact);
bits=bitor(uint32(typecast(reshape(words(1,:),[],1),'uint16')), ...
    bitshift(uint32(typecast(reshape(words(2,:),[],1),'uint16')),16));
iq(double(index)+1)=typecast(bits,'int32');
end
