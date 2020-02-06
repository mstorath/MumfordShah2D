function weights = plExpandWeights( weights, f)
%PLEXPANDWEIGHTS Expand weights for vector-valued data
if size(f,3) ~= size(weights,3)
    weights = repmat(weights, [1 1 size(f,3)]);
end

end

