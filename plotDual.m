% This code accompanies the paper "On the frame set of the second-order B-spline"
% by A. G. D. Atindehou, C. Frederick, Y. B. Kouagou, And K. A. Okoudjou
% https://arxiv.org/abs/1806.05614


function plotDual
a = .85;
b = 1.05;
N = 2;
g = @(x) fnval(spmak(-N:N,[0 1 0]),x);
x=-a/2:.001:a/2;

% Set up and solve the linear systems for h(x)
for j=1:length(x)
    [m, Gm]=G(x(j),g, a, b);
    if m==0
        H = NaN;
    else
        bv = zeros(2*m-1,1); 
        bv(m) = b;
        H(:,j) = Gm\bv;
    end
end

% concatenate the components of H
X = []; h = [];
for k = -(m-1):(m-1)
    X = [X x+k*a];
    h = [h, H(k+m,:)];
end

figure;
plot(X(abs(h)>1e-10), h(abs(h)>1e-10), '.')
title(sprintf('Dual of $\\mathcal{G}(B_N)$ \n N=%d, a=%1.2f, b=%1.2f,',N, a, b))
end

function [m, Y] = G(x, g, a, b, m_max)
% Returns the (2m-1) x (2m-1) matrix Y = G_3(x) and the integer m between 1 and m_max such that (a,b) is in \Lambda_m.
% If (a,b) is not in \Lambda_m for any m\geq 1, Y=NaN.
% The window g is a function handle, m_max is a large integer and x, a, and b are real numbers.

if nargin<5
    m_max=50;
end

m = NaN;
for mm = 1:m_max
    if (b>2*(mm-1)/(2+(2*mm-3)*a)) && (b<2*mm/(2+(2*mm-1)*a)) && (b<2/(1+a)) && (b>1)
         m= mm;
    end
end
if ~isnan(m)
    l = meshgrid(-(m-1):(m-1));
    Y = g(x-l'/b+l*a);
else
    Y=NaN;
    fprintf('Error. Perhaps try again with larger m_max.')
end
end
