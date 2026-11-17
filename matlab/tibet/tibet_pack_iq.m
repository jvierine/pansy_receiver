function [compact,index,words]=tibet_pack_iq(iq)
% Exact int32 storage using int16 plus sparse original low/high words.
assert(isa(iq,'int32'),'Expected original integer voltages');
outside=iq<intmin('int16') | iq>intmax('int16');
compact=int16(iq);
index=uint64(find(outside)-1);
original=iq(outside); bits=typecast(original(:),'uint32');
lo=typecast(uint16(bitand(bits,uint32(65535))),'int16');
hi=typecast(uint16(bitshift(bits,-16)),'int16');
words=[lo(:)';hi(:)'];
assert(isequal(iq,tibet_unpack_iq(compact,index,words)),'IQ round trip failed');
end
