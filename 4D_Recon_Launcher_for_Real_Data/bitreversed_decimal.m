function [out]=bitreversed_decimal(in, maxbits)

if(maxbits == 0)
    out = 0;
    return;
end

in=dec2bin(in, maxbits);
out(1:maxbits)='0';
for i=1:maxbits
    out(i)=in(maxbits-i+1);
end
out=bin2dec(out);