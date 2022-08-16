function [varargout] = dagaf(S, nIMF, chi, eType1, eType2, TH1, TH2)
% Data-adaptive Gaussian average filtering on signal S.
%
% Designed by Yue-Der Lin and Kai-Chun Lin, Taiwan.
% Reference:
%    Yue-Der Lin*, Yong Kok Tan  and Baofeng Tian, "Mode decomposition of 
%    biomedical signals by data-adaptive Gaussian averaging filtering", 
%    Biomedical Signal Processing and Control, Vol.71, Part A, 103104, 
%    January 2022. https://doi.org/10.1016/J.BSPC.2021.103104
%
% Welcome to use this program, and please cite the above paper. 
%
% Inputs:
%       S:    The signal to be decomposed (must be an array).
%    nIMF:    Number of IMF expected to decompose (default = 2).
%     chi:    The parameter that influence the window length,
%             (default = 1.6). The value should be between 1.1 and 3.
%    eType1:  The extention type for the left boundary point:
%             'p' (periodical) repeats the pattern outside the boundaries,
%             'c' (constant) extends outside the boundaries with the last
%                  values achieved at the boundaries,
%             'r' (reflection) extends the signal symmetrical with respect
%                  to the vertical lines over the boundary points.
%             'd' (double-symmetric reflection) extends the signal firstly
%                  symmetrical with respect to the vertical lines over the 
%                  boundary point and and next symmetrical with respect to
%                  horizontal line over the boundary point (default).
%    eType2:  The extention type for the right boundary point:
%             'p' (periodical) repeats the pattern outside the boundaries,
%             'c' (constant) extends outside the boundaries with the last
%                  values achieved at the boundaries,
%             'r' (reflection) extends the signal symmetrical with respect
%                  to the vertical lines over the boundary points.
%             'd' (double-symmetric reflection) extends the signal firstly
%                  symmetrical with respect to the vertical lines over the 
%                  boundary point and and next symmetrical with respect to
%                  horizontal line over the boundary point (default).
%    TH1:     threshold value for signal to residual energy ratio,
%             the 1st decomposition stop criteria (default = 20), and is
%             computed as 10*log10(||S(t)||^2/||S_k(t)||^2), where
%                         S(t) is the original signal, 
%                         S_k(t) is the residual of the kth IMF.
%    TH2:     threshold value for convergence check,
%             the 2nd decomposition stop criteria (default = 0.001), and is
%             computed as ||{S_{i-1}(t)-S_{i}(t)||^2/||S_{i}(t)||^2 at i-th
%             iteration.
% Outputs:
%    imf:     Matrices containg in row j the j-th imf. 
%    res:     S - (summation of imf).
% Called functions:
%    numExtrema.m:  To find the number of extremes in signal S.
% Example:
%    >> load LOD.txt;
%    >> [imf, res] = dagaf(LOD, 4, 1.6, 'd', 'd', 20, 0.001);

% Default arguments:
% Inputs:
if nargin == 1 
    nIMF=2; chi=1.6; eType1='d'; eType2='d'; TH1=20; TH2=0.001;
end
if nargin == 2 
    chi=1.6; eType1='d'; eType2='d'; TH1=20; TH2=0.001;
end
if nargin == 3 
    eType1='d'; eType2='c'; TH1=20; TH2=0.001;
end
if nargin == 4 
    eType2='c'; TH1=20; TH2=0.001;
end
if nargin == 5 
    TH1=20; TH2=0.001;
end
if nargin == 6
    TH2=0.001;
end
if nargin > 7
    disp('Too many arguments.\n');
end
  
% For the cases of wrong inputs:
% For eType1:
if eType1 ~= 'p' || eType1 ~= 'c' || eType1 ~= 'r' || eType1 ~= 'd'
    eType1 = 'd';
end
% For eType2:
if eType2 ~= 'p' || eType2 ~= 'c' || eType2 ~= 'r' || eType2 ~= 'd'
    eType2 = 'd';
end
% For nIMF:
if nIMF <= 1 || nIMF ~= floor(nIMF)
    nIMF = 2;
end
% For chi:
if chi > 3 || chi < 1.1
    chi = 1.6
end
% For TH1 and TH2:
if TH1 < 0 
    TH1 = 1;
end
if TH2 < 0 || TH2 >= 1
    TH2 = 0.001;
end

% Initialization:
L = length(S); S = reshape(S, 1, L); % Make it a row array.
% Initialize IMF to be a null matrix.
imf = [];
    
% energyRatio: pre-defined threshold for signal to residual energy ratio.
energyRatio = 10; % Sifting as energyRatio < TH1:
% rTolerancxe: pre-defined threshold for convergence check
rTolerance = 1;  % Sifting as rTolerance > TH2:

S1 = S;

% DAGAF:
% ind: index for IMF matrix, indicating the ind-th row.
ind = 1; 

% Plot the average filtering process:
figure; 

while (ind <= nIMF) % Outer loop:
    % The first convergence check:
    energyRatio = 10*log10(norm(S,2)/norm(S1,2));
    
    [MaxMin, ~, ~] = numExtrema(S1);
    nMaxMin = length(MaxMin);
    
    if (energyRatio > TH1) || (nMaxMin <= 2)
        break
    end
    
    % Sifting process initialization:
    S2 = S1;
    
    % Sifting process:
    rsigPrev = S2;

    % Filtering procedure:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [MaxMin, Max, Min] = numExtrema(S1);
    [MaxMin, ~, ~] = numExtrema(S1);
    nMaxMin = length(MaxMin);
        
    if nMaxMin > 2
        % You can change the parameter 1.6 to be a number within 1.1 and 3.
        mask = 2*floor(chi*L/nMaxMin);
        H = gausswin(2*mask+1, 4.0728); H = reshape(H, 1, 2*mask+1);
        W = H/sum(H);
        
        % Generate the new pattern according to eType1 and eType2:
        % For left boundary point:
        if eType1 == 'p'
            St_tmp = [S2(end-mask:end-1) S2];
        elseif eType1 == 'r'
            St_tmp = [fliplr(S2(2:mask+1)) S2];
        elseif eType1 == 'd'
            xt1 = 2*mean(S2) - fliplr(S2(2:mask+1));
            St_tmp = [xt1 S2];
        else % eType1 == 'c'
            St_tmp = [S2(1).*ones(1,mask) S2];
        end
            
        % For right boundary point:
        if eType2 == 'p'
            St = [St_tmp S2(2:mask+1)];
        elseif eType2 == 'r'
            St = [St_tmp fliplr(S2(end-mask:end-1))];
        elseif eType2 == 'd'
            xt2 = 2*mean(S2) - fliplr(S2(end-mask:end-1));
            St = [St_tmp xt2];
        else % eType2 == 'c'
            St = [St_tmp S2(end).*ones(1,mask)];
        end
            
        % Filtering:
        for i = 1 : L
            ave(i) = W*St(i:i+2*mask)';
        end
    else
        ave = zeros(1, L);
    end
        
    ave = reshape(ave, 1, L);
    subplot(nIMF,1,ind); plot([1:L],S2,'k',[1:L],ave,'r','LineWidth', 1.5);
    % title('Signal (black) and Average (red)')
    set(gca,'FontSize',10,'FontName','palatino linotype','FontWeight','bold');
    S2 = S2 - ave;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    % Residual tolerance:
    rTolerance = (norm(rsigPrev-S2,2)/norm(rsigPrev,2))^2;
        
    % The second convergence check:
    if (rTolerance < TH2)
        break;
    end
    
    % After sifting:
    imf = [imf; S2];
    S1 = S1 - imf(end,:);
    
    ind = ind + 1;
end

% Plot each inf:
figure; 
[nIMF, ~] = size(imf);
res = S;
for i = 1 : nIMF
    subplot(nIMF, 1, i); plot(imf(i,:),'LineWidth', 1.5); 
    set(gca,'FontSize',10,'FontName','palatino linotype','FontWeight','bold');
    res = res -imf(i,:);
end

% Plot the signal and the residual:
figure; 
plot([1:L],S,'k',[1:L],res,'r','LineWidth', 1.5);
set(gca,'FontSize',14,'FontName','palatino linotype','FontWeight','bold');
% title('Signal (black) and Residual (red)')

% Outputs:
if nargout > 0
    varargout{1} = imf;
end
if nargout > 1
    varargout{2} = res;
end

end % End of dagaf.m

function [MaxMin, Max, Min] = numExtrema(S)
% Find the number of extremas in signal S.
% Inputs:
%    S: the signal to be analyzed (should be an array).
% Outputs:
%    MaxMin: the array for the indices of extrema points.
%    Max:    the array for the indices of maxima points.
%    Min:    the array for the indices of minima points.

% Find the indices of maxima points:
Max    = find(diff(diff(S) > 0) < 0);
ind    = find(S(Max+1) > S(Max));
Max(ind) = Max(ind)+1;

% Find the indices of minima points:
X = -S;
Min    = find(diff(diff(X) > 0) < 0);
ind    = find(X(Min+1) > X(Min));
Min(ind) = Min(ind)+1;

% Find the indices of extrema points:
MaxMin = sort([Max Min]);

end % End of numExtrema.m 