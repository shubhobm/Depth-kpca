function z = kerneldenoise(X,X1,evec,beta,sigma,opt)
% This function is part of the pipeline that produces the results of the
% analysis of the object recognition data set presented in the paper
% entitled "Nonlinear denoising and analysis of 
% neuroimages with kernel principal component analysis and pre-image estimation".
%
% ------------------------------------------------------------------------%
% Trine Abrahamsen, DTU Informatics
% email: tjab@imm.dtu.dk
% Peter Mondrup Rasmussen, DTU Informatics
% email: peter.mondrup@gmail.com  homepage: www.petermondrup.com
%
% ------------------------------------------------------------------------%
% version history:
% Oct 11 2011 - first implementation.
% ------------------------------------------------------------------------%

% Licence: The code is availabel under the MIT License (MIT). See the file
% "licence.txt" for further information.


opt.N = size(X1,2);
z = nan(opt.P,opt.N);
for cv1 = 1:opt.N
    options.print = 0;
    if rem(cv1,5)==0
        fprintf('working on image %d of %d\n',cv1,opt.N);
        options.print = 1;
    end
    %
    %
    switch opt.denoisemethod
        case 'mika'
            if options.print; fprintf('Denoising with Mika''s method\n'),end
            randn('state',1);rand('state',1); % to reproduce results in paper
            options.tmax = 200;
            options.eps = 1e-9;
            init_point = X1(:,cv1); % use observation as initial point
            z(:,cv1) = mikaaux(evec,beta(:,cv1), X,init_point,sigma,options);
            %
        case 'kwok'
            if options.print; fprintf('Denoising with Kwok and Tsang''s method\n'),end
            options.nn = opt.nn;
            kopt.type = 'rbf';
            kopt.sigma = sigma;
            if cv1 == 1
                [K, K1] = kernel(X,X1,kopt);
            end
            k = K1(:,cv1); % select kernel elements corresponding to the image to be denoised
            z(:,cv1) = kwokaux(evec,X,K,k,sigma,options);
    end
end
end
%
%% Define preimage estimation functions
%
function z = mikaaux(evec,beta,X,init_point,sigma,options)
% This function estimates preimage based on Mika's method.
% "Mika, S., Scholkopf, B., Smola, A., Muller, K. R., Scholz, M., Ratsch, G.,
% Kernel PCA and de-noising in feature spaces,
% Advances in Neural Information Processing Systems 11,
% MIT Press, pp. 536–542, 1999."
%
%
% Input:    alpha       -   eigenvectors (N x K)   (scaled, 1 = nu*alpha'*alpha and sorted)
%           beta       -   projection of the point to be denoised. alpha'*Ktest
%           X           -   input point associated with each coeff [P x N]
%           sigma       -   width of the rbf kernel
%           init_point  -   initial point for optimization
%
% ------------------------------------------------------------------------%
% Trine Abrahamsen, DTU Informatics
% email: tjab@imm.dtu.dk
% Peter Mondrup Rasmussen, DTU Informatics
% email: pmra@imm.dtu.dk
%
% ------------------------------------------------------------------------%
% version history:
% Oct 11 2011 - first implementation.
% ------------------------------------------------------------------------%
%
N = size(X,2);
z = init_point;
old_z = z;
iter = 0;
change = inf;
ch = nan;
options.printlevel = 1;
%
% set options for kernel calculation
kopt.type = 'rbf';
kopt.sigma = sigma;
kopt.test = 1; % only calcluate test kernel
%
% calculate coefficients
gamma = evec*beta;
gammatilde = gamma(:) + (1/N)*(1-sum(gamma(:)));
%
% main loop
while iter < options.tmax && change > options.eps,
    iter = iter + 1;
    [~ , k] = kernel(X,z,kopt);
    ka = k.*gammatilde(:);
    den = sum(ka);
    z = X*ka/den;
    change = sum((old_z-z).^2)/sum((old_z).^2);
    ch(iter) = change;
    if options.print && rem(iter,options.printlevel) == 0
        fprintf('Mika iteration %d, change %e\n',iter,change);
    end
    old_z = z;
end
end
%
function z = kwokaux(evec,X,K,k,sigma,options)
% This function estimates preimage based on Kwok and Tsang's method.
% "Kwok, J. T.-Y., Tsang, I. W.-H.
% The pre-image problem in kernel methods.
% IEEE transactions on neural networks 15 (6), 1517–1525, 2004."
%
%
% Input:    evec        -   eigenvectors (N x K)   (scaled, 1 = nu*alpha'*alpha and sorted)
%           X           -   input point associated with each coeff [P x N]
%           K           -   Kernel matrix holding examples (N x N)
%           k           -   Kernel vector holding elements for example to be denoised (N x 1)
%           sigma       -   width of the rbf kernel
%           options.nn  -   number of neighbours
%
% ------------------------------------------------------------------------%
% Trine Abrahamsen, DTU Informatics
% email: tjab@imm.dtu.dk
% Peter Mondrup Rasmussen, DTU Informatics
% email: pmra@imm.dtu.dk
%
% ------------------------------------------------------------------------%
% version history:
% Oct 11 2011 - first implementation.
% ------------------------------------------------------------------------%


%% Some constants
nn = options.nn;
%
N = size(X,2);
H = eye(N)-1/N*ones(N);
M = evec*evec';
o = ones(N,1);
%
%
%% Feature space distances squared
q = H'*M*H*(k-1/N*K*o);
const = (k+1/N*K*o)'*q + 1/power(N,2)*o'*K*o + 1;
df2 = -2*K*q - 2/N*K*o + const; % eq. 9 in Kwok and Tsang's paper
%
%% Input space distance squared
d2 = real(-sigma * log( 1 - 0.5*df2));
%
%% Select nn neighbours
[~, inx] = sort(df2);
X = X(:,inx(1:nn));
d2 = d2(inx(1:nn));
H = eye(nn,nn) - 1/nn * ones(nn,nn);
[U,L,V] = svd(X*H,0);
r = rank(L,1e-5*max(diag(L)));
U = U(:,1:r);
L = L(1:r,1:r);
V = V(:,1:r);
Z = L*V';
d02 = sum(Z.^2)';
z = -0.5*pinv(Z')*(d2-d02);
z = U*z + sum(X,2)/nn;
end