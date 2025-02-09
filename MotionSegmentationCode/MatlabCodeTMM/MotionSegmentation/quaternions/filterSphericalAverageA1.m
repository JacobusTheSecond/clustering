function [Y,max_it,t] = filterSphericalAverageA1(varargin)
% Y = filterSphericalAverageA1(w,X,step,padding_method)
% Uses the function sphericalAverageA1 to filter curves embedded in
% the unit quaternion sphere with a sliding window.
%
% Input:    w, weight vector
%           X, 4xN matrix of input unit quaternions
%           optional: step, step size for window (window length is length(w)), default is step=1.
%               ->!! Length of output sequence is ceil(N/step). !! <-
%           optional: padding method, one of {'symmetric', 'zero'}. default is 'symmetric'
%
% Output:   Y, Filtered version of X
%           max_it, maximum number of iterations for computations of minimum weighted distances
%           t, running time for filter.
%
% Reference: S. Buss, J. Fillmore: Spherical Averages and Applications to
% Spherical Splines and Interpolation, ACM Trans. Graphics 20 (2001):
% 95--126


switch (nargin)
    case 2
        w = varargin{1};
        X = varargin{2};
        step = 1;
        padding_method = 'symmetric';
    case 3
        w = varargin{1};
        X = varargin{2};
        step = varargin{3};
        padding_method = 'symmetric';
    case 4
        w = varargin{1};
        X = varargin{2};
        step = varargin{3};
        padding_method = varargin{4};
    otherwise
        error('Wrong number of arguments!');
end

L = size(w,2);
N = size(X,2);

if (L>N)
    error('Filter length must not be larger than number of data points!');
end

%%%% prepare data set X by means of pre- and postpadding
switch mod(L,2)
    case 0 % even filter length
        prepad_length = L/2;
        postpad_length = L/2 - 1;
    case 1 % odd filter length
        prepad_length = (L - 1)/2;
        postpad_length = (L - 1)/2;
end

tic;
if (strncmp(padding_method,'symmetric',1))
    pre = fliplr(X(:,1:prepad_length));
    post = fliplr(X(:,N-postpad_length+1:N));
elseif (strncmp(padding_method,'zero',1))
    pre = [ones(1,prepad_length);zeros(3,prepad_length)];
    post = [ones(1,postpad_length);zeros(3,postpad_length)];;
else
    error('Unknown padding option!');
end
X = [pre X post];

Y = zeros(4,ceil(N/step));
it = 0; max_it = 0;
for (i=1:ceil(N/step))
    pos = step*(i-1)+1;
    [Y(:,i),it] = sphericalAverageA1(X(:,pos:pos+L-1),w);
    if (it > max_it)
        max_it = it;
    end
end
t = toc;