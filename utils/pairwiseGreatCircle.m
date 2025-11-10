function D = pairwiseGreatCircle(lat, lon)
%PAIRWISEGREATCIRCLE  Vectorised haversine distance matrix (km).
%
%   D = pairwiseGreatCircle(lat, lon) returns an n‑by‑n symmetric matrix of
%   great‑circle distances between points given by (lat,lon) in decimal
%   degrees.  Uses the mean Earth radius R = 6 371 km.
%
%   lat, lon  column vectors of equal length n.

R = 6371;                              % mean Earth radius [km]
lat  = deg2rad(lat);                   % radians
lon  = deg2rad(lon);

% Compute all pairwise differences via implicit expansion
dlat = lat - lat.';                    % n‑by‑n
dlon = lon - lon.';

sin2  = sin(dlat/2).^2 + ...
        cos(lat).*cos(lat.').*sin(dlon/2).^2;

% Numerical safety: asin argument ≤1
sin2 = min(1, max(0, sin2));
c    = 2 * asin( sqrt(sin2) );

D    = R * c;                          % [km]
end