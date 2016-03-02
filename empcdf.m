function [x, p] = empcdf(x)
%EMPCDF Empirical cumulative distribution function.
%
%   [X, P] = EMPCDF(X) is the empirical cumulative distribution function of the
%   sample in X, P(i) = P{X <= x(i)}.  X must be a vector.
%
%   For instance, to plot the empirical cdf for two small equal-size
%   samples, one may use
%
%      x = randn(10, 1);            % generate sample
%      [x, p] = empcdf(x);          % compute empirical cdf
%      stairs(x, p);                % plot the cdf
%      set(gca, 'YLim', [0 1]);     % make sure y-range is [0,1]
%
%   See also ECDF (Statistics Toolbox).

%   Author:      Peter John Acklam
%   Time-stamp:  2004-02-22 12:05:12 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

   % Check number of input arguments.
   nargsin = nargin;
   error(nargchk(1, 2, nargsin));

   sx = size(x);                % size of `x'
   dx = ndims(x);               % number of dimensions in `x'

   % Get first non-singleton dimension, or 1 if none found.
   if nargsin < 2
      k = find(sx ~= 1);
         if isempty(k)
         dim = 1;
      else
         dim = k(1);
      end
   else
      if any(size(dim) ~= 1) | (dim < 1) | (dim ~= fix(dim))
         error('Dimension must be a scalar positive integer.');
      end
   end

   n = size(x, dim);
   x = sort(x, dim);

   % Create a vector P = [1/N, 2/N, ..., 1] along dimension DIM.
   p = (1 : n) / n;
   siz = ones(1, dx);
   siz(dim) = n;
   p = reshape(p, siz);

   % Replicate P along all dimensions except DIM to match the size of X.
   rep = sx;
   rep(dim) = 1;
   p = repmat(p, rep);