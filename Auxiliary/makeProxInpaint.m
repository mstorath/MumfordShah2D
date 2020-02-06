function handle = plMakeProxInpaint( f, mask)

mask = plExpandWeights(mask, f);
handle = @(z, tau, opts) f .* mask + z .* ~mask;
    
end




