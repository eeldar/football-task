classdef Utilities < handle
    
    properties
        
    end
    
    methods(Static)

        function [posp, negp, possumall, posstartall, poslengthall, negsumall, negstartall, neglengthall] = timeseries(series, numclusters, thresh)
            % INPUT:  series - K x J x I data matrix with K subjects, J timepoints, and I conditions (I = 1 for testing against zero, and I = 2 for comparing between conditions) 
            %         numclusters - integer (optional; default = 1)
            %         indicating which supra-threshhold cluster to consider (1 - the one with highest sum of t values, 2 - the second highest sum of t values, etc.)
            %         thresh - scalar (optional; default = 0) indicating t value to use as threshold for identifying clusters
            %
            % OUTPUT: posp - one-tailed p value for condition 1 higher than condition 2 (or greater than zero if there is only one condition)
            %         negp - one-tailed p value for condition 1 lower than condition 2 (or zero)
            %         possumall - sum of t values for the posp cluster 
            %         posstartall - starting timepoint of the posp cluster 
            %         poslengthall - number of timepoints in the posp cluster 
            %         negsumall - sum of t values for the negp cluster 
            %         negstartall - starting timepoint of the negp cluster 
            %         neglengthall - number of timepoints in the negp cluster 
            %
            % EXAMPLE: Utilities.timeseries(randn(10,100,2))
            %
            if nargin<2 || isempty(numclusters); numclusters = 1; end
            if nargin<3; thresh = 0; end
            
            for si =1:size(series,4)
                if size(series,3)==1
                    [~, ~, ~, stats] = ttest(series(:,:,:,si));                
                else
                    [~, ~, ~, stats] = ttest2(series(:,:,1,si),series(:,:,2,si));
                end
                t = stats.tstat;
                [possum, negsum, posstart, poslength, negstart, neglength] = Utilities.consec(t,thresh);
                if length(posstart)>=numclusters
                    possumall(si) = possum(numclusters);
                    possumallsum(si) = sum(possum(1:numclusters));
                    posstartall(si) = posstart(numclusters);
                    poslengthall(si) = poslength(numclusters);
                else
                    possumall(si) =nan;
                    possumallsum(si) = nan;
                    posstartall(si) = nan;
                    poslengthall(si) = nan;
                end
                if length(negstart)>=numclusters
                    negsumall(si) = negsum(numclusters);
                    negsumallsum(si) = sum(negsum(1:numclusters));
                    negstartall(si) = negstart(numclusters);
                    neglengthall(si) = neglength(numclusters);
                else
                    negsumall(si) =nan;
                    negsumallsum(si) =nan;
                    negstartall(si) = nan;
                    neglengthall(si) = nan;
                end
            end
            
            for rep = 1:10000
                if size(series,3)==1
                    ord = repmat(round(rand(size(series,1),1))*2-1,[1 size(series,2)]);
                    for si = 1:size(series,4)
                        [~, ~, ~, stats(si)] = ttest(series(:,:,:,si).*ord);
                    end
                else
                    ord = randperm(size(series,1)*2);
                    while all(ord(1:size(series,1))<=size(series,1)) || all(ord(1:size(series,1))>size(series,1))
                        ord = randperm(size(series,1)*2);
                    end
                    ord1 = ord(ord(1:size(series,1))<=size(series,1));
                    ord2 = ord(ord(1:size(series,1))>size(series,1)) - size(series,1);
                    for si = 1:size(series,4)
                        rseries(:,:,1) = cat(1,series(ord1,:,1,si),series(ord2,:,2,si));                    
                        rseries(:,:,2) = cat(1,series(setdiff(1:size(series,1),ord1),:,1,si),series(setdiff(1:size(series,1),ord2),:,2,si));                    
                        [~, ~, ~, stats(si)] = ttest2(rseries(:,:,1),rseries(:,:,2));
                    end
                end
                for si = 1:size(series,4)
                    t = stats(si).tstat;
                    [possum, negsum, posstart, ~, negstart] = Utilities.consec(t,thresh);
                    if length(posstart)>=numclusters
                        rpossum(si,rep) =possum(numclusters);
                        rpossumsum(si,rep) =nansum(possum(1:numclusters));
                    else
                        rpossum(si,rep) = 0;
                        rpossumsum(si,rep) = nansum(possum);
                    end
                    if length(negstart)>=numclusters
                        rnegsum(si,rep) = negsum(numclusters);
                        rnegsumsum(si,rep) =nansum(negsum(1:numclusters));
                    else
                        rnegsum(si,rep) = 0;
                        rnegsumsum(si,rep) = nansum(negsum);
                    end
                end
            end
            posp = nanmean(any(rpossum>repmat(possumall,[size(series,4) 1]) | rpossumsum>repmat(possumallsum,[size(series,4) 1]),1));
            negp = nanmean(any(rnegsum<repmat(negsumall,[size(series,4) 1]) | rnegsumsum<repmat(negsumallsum,[size(series,4) 1]),1));
            if isnan(posstartall); posp = nan; end
            if isnan(negstartall); negp = nan; end

        end
        
        function [possum, negsum, posstart, poslength, negstart, neglength] = consec(series, thresh)
            if nargin<2; thresh = 0; end
            
            sgn = series<=thresh | isnan(series);
            sgn = [1 sgn 1];
            poslength = find(diff(sgn)==1)-find(diff(sgn)==-1);
            posstart = find(diff(sgn)==-1);
            for i = 1:length(posstart)
                possum(i) = sum(series(posstart(i):posstart(i)+poslength(i)-1));
            end
            try
                [~, I] = sort(possum, 'descend');
                possum = possum(I);
                poslength = poslength(I);
                posstart = posstart(I);
            catch
                possum = [];
                poslength = [];
                posstart = [];
            end
            
            sgn = series>=-thresh | isnan(series);
            sgn = [1 sgn 1];
            neglength = find(diff(sgn)==1)-find(diff(sgn)==-1);
            negstart = find(diff(sgn)==-1);
            for i = 1:length(negstart)
                negsum(i) = sum(series(negstart(i):negstart(i)+neglength(i)-1));
            end
            try
                [~, I] = sort(negsum, 'ascend');
                negsum = negsum(I);
                neglength = neglength(I);
                negstart = negstart(I);
            catch
                negsum = [];
                neglength = [];
                negstart = [];
            end
        end
                           
        function h = lineplot(sig,cols,x)
            % INPUTS: sig - I x J x K data matrix with I conditions, J timepoints, and K subjects
            %         cols - 1 x K vector (optional) indicating color for each condition, between 0 and 1
            %         x - 1 x J vector (optional) indicating x axis values
            %
            % EXAMPLE: Utilities.lineplot(cat(1,rand(1,100,100),rand(1,100,100)+0.1),[0.6 0.8],.1:.1:10)
            %
            hold all;
            N = size(sig,1);
            for i = 1:N
                sigi = nanmean(sig(i,:,:),3);
                n = (length(sigi)-1)/2;
                semi = nanstd(sig(i,:,:),[],3)./sqrt(size(sig,3));
                if nargin<2 || isempty(cols); col = (i-0.9)/N; else col = cols(i); end
                if nargin<3; x = (0:n*2)/100; end;
                x = x(1:length(sigi));
                Utilities.ciplot(sigi-semi, sigi+semi,x,Utilities.hsv2rgb([col, .5 .9]));
            end            
            for i = 1:N
                sigi = nanmean(sig(i,:,:),3);
                n = (length(sigi)-1)/2;
                if nargin<2 || isempty(cols); col = (i-0.9)/N; else col = cols(i); end
                h(i) = plot(x, sigi, 'Color', Utilities.hsv2rgb([col, 0.9 1]), 'LineWidth', 2);
            end
        end
        
        function ciplot(lower,upper,x,colour)

            % ciplot(lower,upper)       
            % ciplot(lower,upper,x)
            % ciplot(lower,upper,x,colour)
            %
            % Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
            % l and u must be vectors of the same length.
            % Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
            % can be overlayed without a problem. Make them transparent for total visibility.
            % x data can be specified, otherwise plots against index values.
            % colour can be specified (eg 'k'). Defaults to blue.

            % Raymond Reynolds 24/11/06

            if length(lower)~=length(upper)
                error('lower and upper vectors must be same length')
            end

            if nargin<4
                colour='b';
            end

            if nargin<3 || isempty(x)
                x=1:length(lower);
            end

            % convert to row vectors so fliplr can work
            if find(size(x)==(max(size(x))))<2
            x=x'; end
            if find(size(lower)==(max(size(lower))))<2
            lower=lower'; end
            if find(size(upper)==(max(size(upper))))<2
            upper=upper'; end

            fill([x fliplr(x)],[upper fliplr(lower)],colour, 'FaceAlpha', 0.4, 'EdgeColor',colour)


        end
        
        function [rout,g,b] = hsv2rgb(hin,s,v)
            %HSV2RGB Convert hue-saturation-value colors to red-green-blue.
            %   M = HSV2RGB(H) converts an HSV color map to an RGB color map.
            %   Each map is a matrix with any number of rows, exactly three columns,
            %   and elements in the interval 0 to 1.  The columns of the input matrix,
            %   H, represent hue, saturation and value, respectively.  The columns of
            %   the resulting output matrix, M, represent intensity of red, blue and
            %   green, respectively.
            %
            %   RGB = HSV2RGB(HSV) converts the HSV image HSV (3-D array) to the
            %   equivalent RGB image RGB (3-D array).
            %
            %   As the hue varies from 0 to 1, the resulting color varies from
            %   red, through yellow, green, cyan, blue and magenta, back to red.
            %   When the saturation is 0, the colors are unsaturated; they are
            %   simply shades of gray.  When the saturation is 1, the colors are
            %   fully saturated; they contain no white component.  As the value
            %   varies from 0 to 1, the brightness increases.
            %
            %   The colormap HSV is hsv2rgb([h s v]) where h is a linear ramp
            %   from 0 to 1 and both s and v are all 1's.
            %
            %   See also RGB2HSV, COLORMAP, RGBPLOT.

            %   Undocumented syntaxes:
            %   [R,G,B] = HSV2RGB(H,S,V) converts the HSV image H,S,V to the
            %   equivalent RGB image R,G,B.
            %
            %   RGB = HSV2RGB(H,S,V) converts the HSV image H,S,V to the 
            %   equivalent RGB image stored in the 3-D array (RGB).
            %
            %   [R,G,B] = HSV2RGB(HSV) converts the HSV image HSV (3-D array) to
            %   the equivalent RGB image R,G,B.

            %   See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78.
            %   Copyright 1984-2011 The MathWorks, Inc. 

            if nargin == 1 % HSV colormap
                threeD = ndims(hin)==3; % Determine if input includes a 3-D array
                if threeD,
                    h = hin(:,:,1); s = hin(:,:,2); v = hin(:,:,3);
                else
                    h = hin(:,1); s = hin(:,2); v = hin(:,3);
                end
            elseif nargin == 3
                if ~isequal(size(hin),size(s),size(v)),
                    error(message('MATLAB:hsv2rgb:InputSizeMismatch'));
                end
                h = hin;
            else
                error(message('MATLAB:hsv2rgb:WrongInputNum'));
            end    

            h = 6.*h;
            k = floor(h);
            p = h-k;
            t = 1-s;
            n = 1-s.*p;
            p = 1-(s.*(1-p));

            % Processing each value of k separately to avoid simultaneously storing
            % many temporary matrices the same size as k in memory
            kc = (k==0 | k==6);
            r = kc;
            g = kc.*p;
            b = kc.*t;

            kc = (k==1);
            r = r + kc.*n;
            g = g + kc;
            b = b + kc.*t;

            kc = (k==2);
            r = r + kc.*t;
            g = g + kc;
            b = b + kc.*p;

            kc = (k==3);
            r = r + kc.*t;
            g = g + kc.*n;
            b = b + kc;

            kc = (k==4);
            r = r + kc.*p;
            g = g + kc.*t;
            b = b + kc;

            kc = (k==5);
            r = r + kc;
            g = g + kc.*t;
            b = b + kc.*n;

            if nargout <= 1
                if nargin == 3 || threeD 
                    rout = cat(3,r,g,b);
                else
                    rout = [r g b];
                end
                rout = bsxfun(@times, v./max(rout(:)), rout);
            else
                f = v./max([max(r(:)); max(g(:)); max(b(:))]);
                rout = f.*r;
                g = f.*g;
                b = f.*b;
            end
        end
        
    end
end
