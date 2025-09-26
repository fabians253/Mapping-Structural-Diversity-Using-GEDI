%% Mapping-Structural-Diversity-Using-GEDI
% 
% main author: Fabian D. Schneider
% 
% This is a function to scale data.
% 
% Please read and reference (cite) the following scientific paper when using this code:
% 
% Fabian D. Schneider*, Morgan Dean, Elsa M. Ordway, Moses B. Libalah, & Antonio A. Ferraz. Mapping the structural diversity of Central African and Western US forests using GEDI. In Review at Remote Sensing of Environment.
% *fabian.schneider@bio.au.dk; Section for Ecoinformatics & Biodiversity, Department of Biology, Aarhus University, Ny Munkegade 114, DK-8000 Aarhus, Denmark

function out = scaleNoOutliers_cols( in, prcInterval )

minPrc = prctile( in, prcInterval );
maxPrc = prctile( in, 100-prcInterval );

for i = 1:size(in,2)
    temp = in(:,i);
    temp( temp < minPrc(i) ) = minPrc(i);
    temp( temp > maxPrc(i) ) = maxPrc(i);
    out(:,i) = temp;
end