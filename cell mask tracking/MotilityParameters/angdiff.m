function d = angdiff(th1, th2)
%%ANGDIFF Difference of two angles
%
% D = ANGDIFF(TH1, TH2) returns the difference between angles TH1 and TH2 on
% the circle.  The result is in the interval [-pi pi).  If TH1 is a column 
% vector, and TH2 a scalar then return a column vector where TH2 is modulo 
% subtracted from the corresponding elements of TH1.
%
% D = ANGDIFF(TH) returns the equivalent angle to TH in the interval [-pi pi).
%
% Return the equivalent angle in the interval [-pi pi).
if nargin < 2
        d = th1;
    else
        d = th1 - th2;
    end
d = mod(d+pi, 2*pi) - pi;
end

