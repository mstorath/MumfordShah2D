function [out,nIter] = mumfordShah2D(gamma, alpha, proxHandle, varargin)
%mumfordShah2D Edge preserving image restoration Mumford-Shah using the
%Mumford-Shah model
% This is an implementation of the method described in the paper
% K. Hohm, M. Storath, A. Weinmann
% An algorithmic framework for Mumford–Shah regularization of inverse problems in imaging
% Inverse Problems 31 (11), 115011

% parse first few options
ip = inputParser;
ip.KeepUnmatched = 1;
addParamValue(ip, 'muSeq', @(k) k^2.01 * 1e-6 ); % by default quadratic + eps progression
addParamValue(ip, 'isotropic', 2); % near isotropic neighborhood
parse(ip, varargin{:});
par = ip.Results;

% set neighborhood
switch par.isotropic
    case 0
        % anisotropic finite differences
        nhood = { [1,0] };
        omega = { 1 };
    case 1
        % finite differences with diagonals
        nhood = { [1,0], [1,1] };
        omega = { sqrt(2) - 1,...
            1 - sqrt(2)/2 };
    case 2
        % finite differences with knight moves
        nhood = { [1,0], [1,1], [2,1], [1,2] };
        omega = { sqrt(5) - 2,...
            sqrt(5) - 3/2 * sqrt(2),...
            0.5 * (1 + sqrt(2) - sqrt(5)),...
            0.5 * (1 + sqrt(2) - sqrt(5)) };
        
    otherwise
        error('Value of isotropic must be 0, 1, or 2')
end
S = 2*numel(nhood);

% parse rest of parameters
%addParamValue(ip, 'nuSeq', @(k) 0.25 * par.muSeq(k) );
addParamValue(ip, 'nuSeq', @(k) par.muSeq(k) * S / nchoosek(S,2) );
addParamValue(ip, 'tol', 1e-3);
addParamValue(ip, 'multiThreading', true);
addParamValue(ip, 'verbose', false);
addParamValue(ip, 'f', 0);
addParamValue(ip, 'maxIter', 50000);
addParamValue(ip, 'maxInnerIter', 500);
addParamValue(ip, 'groundTruth', []);
addParamValue(ip, 'caxis', [0,1]);
addParamValue(ip, 'method', 'L2'); % L1 or L2 variational term
addParamValue(ip, 'numThreads', maxNumCompThreads); 
addParamValue(ip, 'filter', 1);

ip.KeepUnmatched = 0;
parse(ip, varargin{:});
par = ip.Results;

if any(strcmp(ip.UsingDefaults,'nuSeq'))
    par.nuSeq = @(k) par.muSeq(k) * S / nchoosek(S,2);
end

assert(par.tol > 0, 'Stopping tolerance must be > 0.');

% options for iterative linear solvers
linopts.maxit = par.maxInnerIter;

% counts total number of iterations
nIter = 0;

% should rho be computed each iteration
rhoFlag = par.nuSeq(1);

% ATf
backproj = proxHandle(0, par.muSeq(1));

% initialize variables
if(ndims(backproj) == 2)
    backproj = shiftdim(backproj,-1);
    if ~isempty(par.groundTruth)
        par.groundTruth = shiftdim(par.groundTruth,-1);
    end
elseif (ndims(backproj) == 3)
    backproj = shiftdim(backproj,2);
    if ~isempty(par.groundTruth)
        par.groundTruth = shiftdim(par.groundTruth,2);
    end
else
    error('Image must have 2 (gray scale) or 3 (color/multichannel) dimensions')
end

% find out maximum size
[~,m,n] = size(backproj);
maxSize = max(m,n);



u = cell(S,1);
lams = cell(S,1);
if rhoFlag
    rho = cell(S,S);
end
for s = 1:S
    u{s} = backproj;
    lams{s} = zeros(size(backproj));
    if rhoFlag
        for t = s+1:S
            rho{s,t} = zeros(size(backproj));
        end
    end
end

v = backproj;
linopts.init = v;

switch par.method
    case 'L1'
        mums = mumfordShah.L2L1MumfordShahDirectionProcessor();
    case 'L2'
        mums = mumfordShah.L2L2MumfordShahDirectionProcessorFaster();
end

% set the amount of threads that should be used
mums.setNumThreads(par.numThreads);

err = 2*par.tol;

while(err > par.tol && par.maxIter > nIter)
    nIter = nIter + 1;
    mu = par.muSeq(nIter);
    nu = par.nuSeq(nIter);
    
    z = zeros(size(backproj));
    for s = 1:S
        w = mu*v + lams{s};
        if rhoFlag
            for r = 1:s-1
                w = w + nu*u{r} + rho{r,s};
            end
            for t = s+1:S
                w = w + nu*u{t} - rho{s,t};
            end
        end
        w = w/(mu+nu*(S-1));
        
        curGamma = 2/(mu+nu*(S-1))*omega{ceil(s/2)}*gamma;
        curAlpha = 2/(mu+nu*(S-1))*omega{ceil(s/2)}*alpha;
        
        switch par.method
            case {'L2','L2Gauss'}
                gauss = mumfordShah.GaussL2Mum(maxSize,curAlpha);
                mums.setGauss(gauss)
            case 'L2Rod'
                gauss = mumfordShah.GaussTautRod(maxSize,curAlpha);
                mums.setGauss(gauss)
        end
        
        if(mod(s,2) == 0)
            mums.set(w, curGamma, curAlpha, nhood{ceil(s/2)}, 0);
            u{s} = mums.run();
        else
            mums.set(rotate90(w,1), curGamma, curAlpha, nhood{ceil(s/2)}, 0);
            u{s} = rotate90(mums.run(),-1);
        end
        z = z + (u{s}-lams{s}/mu);
    end
    
    z = z/S;
    
    v = proxHandle(shiftdim(z,1), mu * S, linopts);
    
    % init h for next Tikhonov iteration ("warm start")
    linopts.init = v;
    if ndims(v) == 2
        v = shiftdim(v,-1);
    else
        v = shiftdim(v, 2);
    end
    
    for s = 1:S
        lams{s} = lams{s} + mu*(v-u{s});
    end
    if rhoFlag
        for r = 1:S-1
            for t = r+1:S
                rho{r,t} = rho{r,t} + nu*(u{r}-u{t});
            end
        end
    end
    
    err = norm(u{1}(:) - u{2}(:),2)/ ( norm(u{1}(:),2) + norm(u{2}(:),2));
    
    if par.verbose 
        if (nIter <= 20)
            modulo = 1;
        else
            modulo = 1;
        end
        
        if(mod(nIter,modulo) == 0)
            out = zeros(size(backproj));
            for s = 1:S
                out = out + u{s};
            end
            out = out/S;
            
            fprintf('%i. iteration\n',nIter)
            fprintf('Mu: %f, ', mu)
            fprintf('Average discrepancy: %f', err)
            if ~isempty(par.groundTruth)
                fprintf(', PSNR: %f', plpsnr(par.groundTruth, out));
            end
            fprintf('\n\n')
            figure(999)
            for ind = 1:numel(nhood)
                subplot(2,numel(nhood),ind)
                imshow(shiftdim(u{2*ind-1},1),par.caxis)
                title(['Index: ' num2str(2*ind-1) ' Iteration: ' num2str(nIter)])
                subplot(2,numel(nhood),numel(nhood)+ind)
                imshow(shiftdim(u{2*ind},1),par.caxis)
                title(['Index: ' num2str(2*ind) ' Iteration: ' num2str(nIter)])
            end
            drawnow
        end
    end
end

close(figure(999))

out = zeros(size(backproj));
for s = 1:S
    out = out + u{s};
end
out = shiftdim(out/S,1);