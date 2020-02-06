function handle = makeProxL2Linop( f, A)

handle = @(v, tau, opts) minL2Tikhonov(f, tau, A, 'u0', v);
    
end




