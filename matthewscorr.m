% R = matthewscorr(X, Y)
%
% computes the 'Matthews correlation coefficient' as described on Wikipedia
% (http://en.wikipedia.org/wiki/Matthews_correlation_coefficient)
%
% It's supposedly a special correlation coefficient for binary data (e.g.
% in binary classification) and computes the correlation based on true and
% false postivies and negatives, but the correlation coefficients seem to
% be the same as those computed by Matlab's corrcoef.
%
% in:
%       X   -   binary data set 1
%               [N, nsets] = size
%       Y   -   binary data set 2
%               [N, nsets] = size
% out:
%       R   -   Matthews correlation coefficient between columns of X and Y
%               [1, nsets] = size
% author:
%       Sebastian Bitzer (bitzer@cbs.mpg.de)
function MCC = matthewscorr(template, attempt)

TN = size(intersect(find(template==0),find(attempt==0)),2);
TP = size(intersect(find(template==1),find(attempt==1)),2);
FP = size(intersect(find(template==0),find(attempt==1)),2);
FN = size(intersect(find(template==1),find(attempt==0)),2);

N = TN + TP + FN + FP;

S = (TP + FN) / N;

P = (TP + FP) / N;

% MCC = ((TP / N) - (S * P))/sqrt((P*S)*(1-S)*(1-P));
MCC = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));

% TP = sum(template == 1 & attempt == 1);
% FP = sum(template == 0 & attempt == 1);
% TN = sum(template == 0 & attempt == 0);
% FN = size(template, 1) - TP - FP - TN;
% 
% R = (TP .* TN - FP .* FN) ./ ...
%     sqrt( (TP + FP) .* (TP + FN) .* (TN + FP) .* (TN + FN) );