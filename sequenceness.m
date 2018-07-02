function seq = sequenceness(X, wind, maxgap)
    % 
    %  compute sequenceness between pairs of time series
    %
    %  USAGE: seq = sequenceness(X, wind, maxgap)
    %
    %  INPUTS: 
    %    X - [T x Q x S] data matrix containing S time series for T trials and Q timepoints per trial
    %    wind - length of window to use for each calculation of sequenceness 
    %    maxgap - maximal time lag between time series to consider
    %
    %  OUTPUTS:
    %    seq - [maxgap x T x P x Q-wind] sequencesness for each time lag upto max gap, for each trial, for each pair of time series, for each starting timepoints
    %
    %  Eran Eldar, July 2018
    
    S = size(X,3); 
    pairs = nchoosek(1:S,2);
    P = size(pairs,1);
    for t = 1:size(X,1)                    
        for p =1:P
            Xt = squeeze(X(t,:,:));
            Xtp = [Xt(:,pairs(p,1)),Xt(:,pairs(p,2))];
            for gap = 1:maxgap
                tm = 1:size(X,2)-wind+1-gap;
                tm1 = repmat(tm,[wind, 1]) + repmat((0:wind-1)',[1, size(tm,2)]);
                tm2 = tm1 + gap;
                X1 = zscore(reshape(Xtp(tm1,:), [size(tm1), 2]));
                X2 = zscore(reshape(Xtp(tm2,:), [size(tm2), 2]));
                temp = mean(X1(:,:,1).*X2(:,:,2) - X1(:,:,2).*X2(:,:,1));
                seq(gap,t,p,:) = cat(4,permute(temp,[4 3 1 2]), nan(1,1,1,size(X,2)-wind-size(temp,2)));
            end
        end
    end
end

