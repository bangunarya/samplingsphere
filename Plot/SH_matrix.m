%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function to generate spherical harmonics matrix
%
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s    

function [ A ] = SphHarm( ang, B )
%SphHarm calculates the spherical harmonics for given theta, phi and order

%B
%   uses Condon-Shortley phase convention to calculate the negative orders
%   https://en.wikipedia.org/wiki/Spherical_harmonics
%   http://scipp.ucsc.edu/~haber/ph116C/SphericalHarmonics_12.pdf
theta = ang(:,1);
phi=ang(:,2);
s = size(theta, 1);
numHarm = (B)^2;

A = zeros(numHarm, s);

i = 0;

for l = 0:B-1
    m = (0:l)';
    
    P = legendre(l, cos(theta));
    
    norm = sqrt( (2*l+1)*factorial(l-m) ./ (4*pi*factorial(l+m)) );
    
    N = norm * ones(1, s);
    
    Exp = exp(1i*m*phi');
    
    Y_pos = N .* P .* Exp;
    
    if l ~= 0
        Condon_Shortley = (-1).^m(end:-1:2) * ones(1, s);
        Y_neg = Condon_Shortley .* conj(Y_pos(end:-1:2, :));
    else
        Y_neg = [];
    end
    
    Y = [Y_neg; Y_pos];
    
    A(i+1:i+(2*l+1), :) = Y;
    i = i+(2*l+1);
end

A = A.';

end
