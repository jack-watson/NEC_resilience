function [D, Z] = simulateNodeDamage(G, params, seed)
%SIMULATENODEDAMAGE  Simulate hourly hurricane‑induced node damage.
%
%   D = SIMULATENODEDAMAGE(G, PARAMS) draws a realisation of the
%   zero‑inflated Beta spatio‑temporal damage model described in the
%   accompanying documentation.  The output is an |V|‑by‑T matrix whose
%   (i,t) entry is the fractional damage d_i(t) sustained by station i
%   during hour t.
%
%   D = SIMULATENODEDAMAGE(G, PARAMS, SEED) additionally sets the random
%   number generator to `rng(SEED)` for reproducibility.
%
% ------------------------------------------------------------------------
% INPUTS
%   G       Undirected MATLAB graph.  G.Nodes must contain:
%             * Name : cell array of char station names  (not used here)
%             * lat  : latitude  [decimal degrees]
%             * lon  : longitude [decimal degrees]
%
%   PARAMS  1‑by‑9 vector of scalar model parameters in the **order**
%             [ phi,          ... %  1 AR‑1 temporal persistence  (‑)
%               sigma2Z,      ... %  2 Latent‑field variance      (‑)
%               rho,          ... %  3 Spatial range (km)         (km)
%               alpha0,       ... %  4 Bernoulli‑logit intercept  (‑)
%               alpha1,       ... %  5 Bernoulli‑logit slope      (‑)
%               beta0,        ... %  6 Beta‑mean‑logit intercept  (‑)
%               beta1,        ... %  7 Beta‑mean‑logit slope      (‑)
%               kappa,        ... %  8 Beta concentration (>0)    (‑)
%               T ];              %  9 Number of hourly steps     (h)
%
%   SEED    (optional) scalar integer random‑seed.
%
% OUTPUT
%   D       |V|‑by‑T matrix of damage fractions d_i(t).
%
% ------------------------------------------------------------------------
% The model:
%   * Latent field  Z_t = phi Z_{t-1} + eps_t,    eps_t ~ N(0, sigma2Z Σ)
%   * Σ_{ij}        = exp( -dist_ij / rho )
%   * Occurrence    A_i(t) ~ Bern( logistic(alpha0 + alpha1 Z_i(t)) )
%   * Magnitude     M_i(t) ~ Beta( mu_i(t)*kappa, (1-mu_i(t))*kappa )
%                    where  mu_i(t) = logistic(beta0 + beta1 Z_i(t))
%   * Damage        d_i(t) = A_i(t) * M_i(t)
%
% The function returns the matrix D without computing subsequent component
% integrities y_i(t) (the caller can accumulate damage if desired).
%
% ------------------------------------------------------------------------
% AUTHOR:  <Your Name>, 2025‑06‑12
% ------------------------------------------------------------------------

%% ---------- 0.  RNG state ------------------------------------------------
if nargin > 2
    rng(seed);                         % user‑supplied seed
end

%% ---------- 1.  Unpack problem data -------------------------------------
n  = numnodes(G);                      % number of stations
lat = G.Nodes.lat(:);                  % column vectors
lon = G.Nodes.lon(:);

% Parameter vector (see header for definition)
phi     = params(1);
sigma2Z = params(2);
rho     = params(3);
alpha0  = params(4);
alpha1  = params(5);
beta0   = params(6);
beta1   = params(7);
kappa   = params(8);
T       = params(9);                   % number of hours to simulate

%% ---------- 2.  Pre‑compute distance & covariance ------------------------
% Pairwise great‑circle distances [km]
Dmat = pairwiseGreatCircle(lat, lon);         % n‑by‑n

% Exponential covariance Σ = exp( -d / rho )
Sigma = exp( -Dmat ./ rho );

% Cholesky factor L such that Σ = L*L'
%   (If Σ is nearly singular add small jitter to the diagonal.)
[L,p] = chol(Sigma, 'lower');
if p > 0                     % not positive‑definite → regularise
    diagJitter  = 1e-10;     % minimal perturbation
    Sigma       = Sigma + diagJitter*eye(n);
    L           = chol(Sigma, 'lower');
end

%% ---------- 3.  Allocate output & initialise latent field ---------------
D = zeros(n, T);                       % preallocate damage matrix

% Stationary variance of AR(1):  var(Z) = sigma2Z / (1 - phi^2)
stationStd = sqrt( sigma2Z / max(1e-12, (1 - phi^2)) );

Zprev = stationStd * ( L * randn(n,1) );   % initial Z_0 ~ N(0, Σ*stationVar)

% Numerically stable logistic/inverse-logit function
logistic = @(x) 1 ./ (1 + exp(-x));

%% ---------- 4.  Main simulation loop ------------------------------------
ZArr = zeros(n, T);
for t = 1:T
    % (a) Draw latent field increment eps_t  ~ N(0, sigma2Z Σ)
    eps_t = sqrt(sigma2Z) * ( L * randn(n,1) );
    
    % (b) Update latent field (AR‑1)
    Zt = phi * Zprev + eps_t;
    
    % (c) Compute occurrence probability & draw Bernoulli
    pOcc   = logistic( alpha0 + alpha1 * Zt );   % n‑by‑1
    A      = rand(n,1) < pOcc;                  % 1=damage event
    
    % (d) Compute Beta mean & draw magnitudes
    mu     = logistic( beta0 + beta1 * Zt );     % n‑by‑1
    aBeta  = mu * kappa;                        % shape‑α
    bBeta  = (1 - mu) * kappa;                  % shape‑β
    
    M = zeros(n,1);            % initialise magnitudes
    idx = A;                   % logical index for damaged nodes
    % betarnd requires both shape parameters >0; logistic guarantees 0<mu<1
    if any(idx)
        M(idx) = betarnd(aBeta(idx), bBeta(idx));
    end
    
    % (e) Hourly damage fraction
    D(:,t) = M;                % because M=0 where A=0
    
    % (f) Prepare for next hour
    Zprev = Zt;
    ZArr(:,t) = Zt;
end

Z = ZArr;

end % ---------- end main function ----------------------------------------
