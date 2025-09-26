%% Mapping-Structural-Diversity-Using-GEDI
% 
% main author: Fabian D. Schneider
% 
% This is the main Matlab script for calculating richness, divergence and evenness based on probability density.
% 
% Please read and reference (cite) the following scientific paper when using this code:
% 
% Fabian D. Schneider*, Morgan Dean, Elsa M. Ordway, Moses B. Libalah, & Antonio A. Ferraz. Mapping the structural diversity of Central African and Western US forests using GEDI. In Review at Remote Sensing of Environment.
% *fabian.schneider@bio.au.dk; Section for Ecoinformatics & Biodiversity, Department of Biology, Aarhus University, Ny Munkegade 114, DK-8000 Aarhus, Denmark

function [fric, feve, fdiv] = getFRicFEveFDiv_PDFadapt( data, xi, probThr, minPoints )

% remove NaNs
uniqueSamples = data( sum( isnan(data), 2 ) == 0, : );
nrSamples = size( uniqueSamples, 1 );

if nrSamples > minPoints

    % add flexible trait scaling

    % zoom to distribution and save trait ranges in all dimensions
    traitMins = prctile( uniqueSamples, 0.01 );
    traitMaxs = prctile( uniqueSamples, 99.99 );
    traitRanges = traitMaxs - traitMins;
    traitMins = traitMins - 0.1*traitRanges;
    traitMaxs = traitMaxs + 0.1*traitRanges;
    traitMins( traitMins < 0 ) = 0;
    traitMaxs( traitMaxs > 1 ) = 1;
    traitRanges = traitMaxs - traitMins;

    % rescale distributions locally
    temp = uniqueSamples - repmat( traitMins, size( uniqueSamples, 1 ), 1 );
    uniqueSamples = temp ./ repmat( traitRanges, size( uniqueSamples, 1 ), 1 );

    % save scale factor of total volume
    scaleFac = prod(traitRanges);

    % tre AKDE function with adaptive kernel

    % optimize speed for large data > 300 samples
    if nrSamples > 300
        gam = ceil( nrSamples^(1/2) * 0.5 );
    else
        gam = ceil( nrSamples^(1/2) );
    end
    try
        trait_pdf = akde( uniqueSamples, xi, gam );
        %figure; plot( trait_pdf );

        % define density threshold for richness as percentage of maximum density
        %densityThreshold = quantile( trait_pdf, probThr );

        % threshold
        ind = trait_pdf > probThr; %max(densityGrid) * 0.5;

        % functional richness as percentage of filled trait space, using density threshold
        fric_temp = sum( ind ) / length( ind );
        fric = fric_temp * scaleFac; % scale back to original volume

        % functional evenness, as deviation from a perfectly even density
        trait_pdf_norm = trait_pdf(ind) / sum(trait_pdf(ind));
        eve_ref = repmat( 1/length( trait_pdf_norm ), 1, length( trait_pdf_norm ) );
        feve = sum( min( cat( 1, eve_ref, trait_pdf_norm(:)' ) ) );

        % functional divergence as the abundance weighted deviation from the
        % mean distance from the center of gravity
        coords = xi(ind,:);
        g = mean( coords, 1 );
        sumDistSq = (coords(:,1)-g(:,1)).^2;
        for j = 2:size( xi,2 )
            sumDistSq = sumDistSq + (coords(:,j)-g(:,j)).^2;
        end
        dg = sqrt( sumDistSq );
        d = sum( trait_pdf_norm .* (dg - mean(dg)) );
        d_abs = sum( trait_pdf_norm .* abs(dg - mean(dg)) );
        fdiv = ( d + mean(dg) ) / ( d_abs + mean(dg) );
    catch
        disp(['Warning: Problem in akde in ' num2str(nrDims) ' dimensions and ' num2str(nrSamples) ' samples.']);
        fric = NaN;
        feve = NaN;
        fdiv = NaN;
    end
else
    fric = NaN;
    feve = NaN;
    fdiv = NaN;
end
