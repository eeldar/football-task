function [lik, latents] = lik_football(P,data)
    %   
    %   Likelihood function for a computational model of the football task 
    %   
    %   USAGE: [lik, latents] = lik_football(P,data)
    %
    %   INPUTS:
    %     P - structure of S parameter samples, with the following fields:
    %         .ba1 - [S x 1] player type 1's scoring
    %         .ba2 - [S x 1] player type 2's scoring
    %         .ba3 - [S x 1] player type 3's scoring
    %         .ba4 - [S x 1] player type 4's scoring
    %         .valtype1 - [S x 1] mixes rounded and unrounded scoring for player type 1
    %         .valtype2 - [S x 1] mixes rounded and unrounded scoring for player type 2
    %         .valtype3 - [S x 1] mixes rounded and unrounded scoring for player type 3
    %         .valtype4 - [S x 1] mixes rounded and unrounded scoring for player type 4
    %         .varnocount - [S x 1] variance per minimially processed player 
    %         .timehitnocountvar - [S x 1] effect of trial time on 'varnocount'
    %         .varcounthit - [S x 1] decreases variance for optimally processed players
    %         .w12 - [S x 1] priority weight for second ranked player type
    %         .w13 - [S x 1] priority weight for third ranked player type
    %         .w14 - [S x 1] priority weight for fourth ranked player type
    %         .w4 - [S x 1] player type 4's priority;
    %         .sequence - [S x 3] mixing coefficients for type-based, numerosity-based, and screen-location-based prioritization
    %         .thresh - [S x 1] resources required for optimal processing
    %         .timehitthresh - [S x 1] effect of trial time on '.thresh'
    %         .inattmean - [S x 1] default scoring coefficient 
    %         .counthitcapall - [S x 1] number of players beyond which default scoring is assumed
    %     data - struture of experimental data with the following fields:
    %         .goals - [1 x T] the participant's answer (between 0 and 10 goals)
    %         .stim - [1 x 4 x T] number of players of each type
    %         .dectime - [1 x T] time avalailable for deliberation (1 or 2 seconds)
    %         .distance - [1 x 4 x T] average distacne from center of screen for each player type
    %
    % OUTPUTS:
    %   lik - [S x 1] log-likelihoods
    %   latents - a structure with the following fields:
    %       .goals - [S x T] model's random answer
    %       .goals_max - [S x T] model's most likely answer
    %       .processing_weights - [S x 4 x T] model's resource allocation
    %
    % Eran Eldar, July 2018
    
    [res, dec, decmax, processing_weights] = prob(P, data.goals, data.stim, data.dectime, data.distance);
    lik = nansum(log(res),2);

    latents.goals = dec;
    latents.goals_max = decmax;
    latents.processing_weights = processing_weights;

end

                
function [res, dec, decmax, processing_weights] = prob(P, goals, stim, dectime, distance)

    [comp, processing_weights] = computation(P, stim, dectime, distance);

    var = 0;
    var = var + bsxfun(@times, P.varnocount, comp(:,:,4)); 
    var(:, dectime == 2) = bsxfun(@times, var(:, dectime == 2), P.timehitnocountvar); 
    varcount = bsxfun(@times, P.varnocount .* P.varcounthit, comp(:,:,3)); 
    var = var + varcount;
    comp = sum(comp(:,:,1:2),3);

    if sum(goals>0)>0; res(:,goals>0) = normcdf(repmat(goals(goals>0)+0.5, [size(comp,1), 1]), comp(:,goals>0,1), sqrt(var(:,goals>0,:))) - normcdf(repmat(goals(goals>0)-0.5, [size(comp,1), 1]), comp(:,goals>0,1), sqrt(var(:,goals>0,:))); end
    if sum(goals==0)>0; res(:,goals==0) = normcdf(0.5*ones(size(comp,1),sum(goals==0)), comp(:,goals==0,1), sqrt(var(:,goals==0,:))); end
    res(:,~(goals>=0)) = nan;

    dec = min(max(round(randn(size(comp(:,:,1))).*sqrt(var) + comp(:,:,1)), 0),10);
    decmax = round(min(max(comp(:,:,1), 0),10));

end

function [comp, counted] = computation(P, stim, dectime, distance)
    
    v = attention(P, stim, dectime, distance);

    ba1 = P.ba1; ba2 = P.ba2; ba3 = P.ba3; ba4 = P.ba4;
    valsa = cat(2, ba1*2/1.6, ba2*1/1.2, ba3*1/0.6, ba4*0/0.4);
    valsb = cat(2, ba1*3/3.2, ba2*2/2.4, ba3*1/1.2, ba4*1/0.8);
    valsc = cat(2, ba1*5/4.8, ba2*4/3.6, ba3*2/1.8, ba4*1/1.2);
    vals = repmat(valsa, [1, 1, size(stim,3)]); 
    for s=1:4; vals(:,s,squeeze(stim(1,s,:))==2) = repmat(valsb(:,s), [1 1 sum(stim(1,s,:)==2,3)]); end
    for s=1:4; vals(:,s,squeeze(stim(1,s,:))==3) = repmat(valsc(:,s), [1 1 sum(stim(1,s,:)==3,3)]); end                  

    vals2 = cat(2, ba1, ba2, ba3, ba4);
    vals2 = repmat(vals2, [1, 1, size(stim,3)]); 
    valtype = cat(2,P.valtype1,P.valtype2,P.valtype3,P.valtype4);
    vals = repmat(valtype, [1 1 size(vals,3)]) .* vals + repmat(1-valtype, [1 1 size(vals,3)]) .* vals2;
    s = repmat(stim, [size(v,1), 1, 1]);

    nocountvals = repmat(P.inattmean, [1, size(vals,2), size(vals,3)]); 

    cap = permute(P.counthitcapall, [1 3 2]);
    cap = repmat(cap, [1, 1, size(v,3)/size(cap,3)]);

    ns = (1-v) .* s;
    cs = s;
    postcap = repmat(max(sum(ns,2)-cap,0),[1 size(v,2)]).*ns./repmat(sum(ns,2),[1 size(v,2)]);
    postcap(isnan(postcap)) = 0;
    notcounted = nocountvals .* postcap + vals .* (ns - postcap);
    counted = vals .* v.* cs;

    comp(:,:,1) = squeeze(sum(counted, 2)); % counted values
    comp(:,:,2) = squeeze(sum(notcounted, 2)); % non-counted values
    comp(:,:,3) = squeeze(sum(v.*s, 2)); % counted figures
    comp(:,:,4) = squeeze(sum((1-v) .* s, 2)); % non-counted figures
    comp(:,:,5) = max(squeeze(sum(v.*(s>0), 2)) - 1, 0); % counted figure types (above one)

    if nargout>1; counted = squeeze(v); end
end
                                    
function v = attention(P, stim, dectime, distance)

    w = [ones(size(P.w12,1), 1, size(P.w12,3)) permute(P.w12, [1 3 2]) permute(P.w13, [1 3 2]) permute(P.w14, [1 3 2])];
    w = [w(:,1) w(:,1).*w(:,2) w(:,1).*w(:,2).*w(:,3) w(:,1).*w(:,2).*w(:,3).*w(:,4)];
    w = w./repmat(sum(w,2), [1, 4, 1]);

    I = ones(1, size(stim,2), size(stim,3));
    s = stim;
    for ns = 1:4
        [~,idx]=sort(s,2,'descend');
        firstOne = repmat(idx(:,1,:), [1, size(idx,2), 1]) == repmat(1:4, [size(s,1), 1, size(s,3)]); 
        I(firstOne) = ns;
        s(firstOne) = -s(firstOne) -1;
        if ns>1
            sameind1 = find(firstOne);
            sameind2 = find(firstOneOld);
            sameind = s(sameind1) == s(sameind2);
            I(sameind1(sameind)) = I(sameind2(sameind));
        end
        firstOneOld = firstOne;
    end
    ind1 = repmat((1:size(w,1))',[size(I,2)*size(I,3),1]);
    ind2 = reshape(repmat(I(:)',[size(w,1), 1]), [length(ind1), 1]);
    if size(w,3)>1
        ind3 = ones(size(w,1)*size(I,2),1)*(1:size(w,3)); ind3 = ind3(:);
        wmajority = w(sub2ind(size(w),ind1, ind2, ind3));
    else
        wmajority = w(sub2ind(size(w),ind1, ind2));
    end
    wmajority = reshape(wmajority, size(w,1),size(I,2),size(I,3));      

    psequence = repmat(permute(P.sequence,[1 4 3 2]), [1,size(stim,2),size(stim,3)/size(P.sequence,3)]);
    I = ones(size(distance));
    for ns = 1:4
        I(distance==repmat(min(distance,[],2), [1, size(stim,2),1]) & ~isinf(distance)) = ns;
        distance(distance==repmat(min(distance,[],2), [1, size(stim,2),1]) & ~isinf(distance)) = inf;
    end
    I(stim==0) = 1;
    ind1 = repmat((1:size(w,1))',[size(I(:),1),1]);
    ind2 = reshape((I(:)*ones(1,size(w,1)))',size(w,1)*size(I(:),1),1);

    if size(w,3)>1
        ind3 = ones(size(w,1)*size(I,2),1)*(1:size(w,3)); ind3 = ind3(:); 
        wdistance = w(sub2ind(size(w),ind1, ind2, ind3));
    else 
        wdistance = w(sub2ind(size(w),ind1, ind2));
    end                    
    wdistance = reshape(wdistance, size(w,1),size(I,2),size(I,3));

    rnk = repmat(3:-1:1, [size(w,1), 1, size(w,3)]);
    w4 = permute(P.w4, [1 3 2]);
    rnk4 = repmat(ceil(w4*4), [1, 3, 1]);
    rnk(rnk>= rnk4) = rnk(rnk>= rnk4) + 1;
    rnk(:,4,:) = rnk4(:,1,:);

    I = ones(size(rnk,1), size(rnk,2), size(stim,3));
    s = repmat(rnk, [1, 1, size(stim,3)/size(rnk,3)]);
    s(repmat(stim,[size(s,1), 1, 1])==0) = 0;
    for ns = 1:4
        [~,idx]=sort(s,2,'descend');
        firstOne = repmat(idx(:,1,:), [1, size(idx,2), 1]) == repmat(1:size(stim,2), [size(s,1), 1, size(s,3)]); 
        I(firstOne) = ns;
        s(firstOne) = -1;
    end
    ind1 = repmat((1:size(w,1))',[size(I,2)*size(I,3),1]);
    ind2 = I(:);
    if size(w,3)>1
        ind3 = ones(size(w,1)*size(I,2),1)*(1:size(w,3)); ind3 = ind3(:);
        wsorted = w(sub2ind(size(w),ind1, ind2, ind3));
    else
        wsorted = w(sub2ind(size(w),ind1, ind2));
    end
    wsorted = reshape(wsorted, size(w,1),size(I,2),size(I,3));  

    w = psequence(:,:,:,1).*wsorted + psequence(:,:,:,2).*wmajority + psequence(:,:,:,3).*wdistance;

    s = repmat(stim, [size(w,1), 1, 1]);
    w(s==0) = 0;

    w = w./repmat(sum(w,2),[1,size(w,2),1]);                

    thresh = permute(P.thresh, [1 3 2]);
    thresh = repmat(thresh, [1, size(w,2), size(w,3)/size(thresh,3)]);
    timehitthresh = permute(P.timehitthresh, [1 3 2]);
    timehitthresh = repmat(timehitthresh, [1, size(thresh,2), size(thresh,3)/size(timehitthresh,3)]);
    thresh(:, : , dectime==2) = thresh(:, : , dectime==2) .* timehitthresh(:, : , dectime==2); 

    v = min(w./thresh, 1);
    v = mean(v,4);
    v(s==0) = 0;
end
     
