function handle = makeProxL2w( f, weights)

weights = expandWeights(weights, f);
handle = @(v, tau, opts) (weights .* f + tau * v) ./ (weights + tau);

end




