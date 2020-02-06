function handle = makeProxL0w( f, weights)

weights = expandWeights(weights, f);
handle = @(z, tau, opts) hardThreshold( z - f, sqrt(2 * weights/tau) ) + f;
    
end




