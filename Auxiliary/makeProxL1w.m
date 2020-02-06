function handle = makeProxL1w( f, weights)

weights = plExpandWeights(weights, f);
handle = @(z, tau, opts) plSoftThreshold( z - f, weights/tau ) + f;
    
end




