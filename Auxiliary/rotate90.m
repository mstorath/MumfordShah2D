function rotated = rotate90(f,k)

if(ndims(f) ~= 3)
    error('Wrong number of dimensions')
end

for i = size(f,1):-1:1
    rotated(i,:,:) = rot90(squeeze(f(i,:,:)),k);
end
