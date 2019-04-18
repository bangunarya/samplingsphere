%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Azimuth search function.
% This code is used to obtain azimuth sequence that could achieve the lower
% bound. It uses pattern search or another heuristic optimization provided
% by MATLAB
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 17.04.2019 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s

function [y_ps,Coh_val_ps] = azimuth_search(initial,x,lm )

    %%%%%%%%%%%%%% OPT PROBLEM
    fun = @(x) simple_objective(x,initial,lm);
   
%     %%%%%%%%%%%%%%% MultiStart
%     gs = MultiStart;
%     problem = createOptimProblem('fmincon','x0',x,...
%             'objective',fun);
%     [y_ps,Coh_val_ps] = run(gs,problem,10) ;
%     y_ps=mod(y_ps,2*pi);
%     %%%%%%%%%%%%%%% PatternSearch
   
       [y_ps,Coh_val_ps] = patternsearch(fun,x) ;
       y_ps=mod(y_ps,2*pi);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function to generate spherical harmonics matrix
%
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s    
    
function [Y,As_Leg]=SH_matrix(ang,lm,basistype)

As_Leg=zeros(size(ang,1),size(lm,1));
Y=zeros(size(ang,1),size(lm,1));
for jj=1:size(lm,1);
    order=lm(jj,:);
 
    
    if isequal(basistype,'complex'); %%
        normalisation=sqrt(((2*order(1)+1)/(4*pi))*(factorial(order(1)-order(2))/factorial(order(1)+order(2))));
        spher_harm=AsslegendreP(order(1),order(2),cos(ang(:,1))).*exp(1i*order(2)*ang(:,2));
        Y(:,jj)=normalisation*spher_harm;
        As_Leg(:,jj)=normalisation*AsslegendreP(order(1),order(2),cos(ang(:,1)));
  

    elseif isequal(basistype,'real');
        if order(2)> 0
            exp_val=sqrt(2)*cos(order(2)*ang(:,2));
        elseif order(2) < 0
            exp_val=sqrt(2)*sin(abs(order(2))*ang(:,2));
        else
            exp_val=ones(size(ang,1),1);
        end%
        normalisation=sqrt(((2*order(1)+1)/(4*pi))*(factorial(order(1)-abs(order(2)))/factorial(order(1)+abs(order(2)))));
    
        spher_harm=normalisation.*AsslegendreP(order(1),abs(order(2)),cos(ang(:,1)));
        As_Leg(:,jj) =spher_harm;
        Y(:,jj)=exp_val.*spher_harm;
    else
        error('What is the basis type? (Complex or Real)');
    end
end
end
% Y_conj=conj(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function to generate associated Legendre Polynomials
%
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s    
    
%% Associated Legendre
function [p]=AsslegendreP(l,m,x)

% Calculate the Associated Legendre polynomials. 
if(nargin==2), 
    x=m;
    clear m;
    m=0;
end;


%Some basic error checking of parameters

% Definition states that if |m|>l, the polynomial is 0; so set to 0 rather
% than return an error. Some algorithms depend on this behaviour. Note that
% the original definition requires 0<=m<=l
if(abs(m)>abs(l)), p=zeros(size(x)); err=0; return; end;

% Definition for l<0
if(l<0), l=-l-1; end;

% For m<0, polynomials are proportional to those with m>0, cfnm is the
% proportionality coefficient
cfnm=1;
if(m<0),
    m=-m;
    cfnm=(-1)^m*factorial(l-m)/factorial(l+m);
end;



% Calculate coef of maximum degree in x from the explicit analytical
% formula
cl=(-1)^m*cfnm*factorial(2*l)/((2^l)*factorial(l)*factorial(l-m));
maxcf=abs(cl);
% fprintf('Coef is %.16f\n',cl);
px=l-m;

% Power of x changes from one term to the next by 2. Also needed for
% sqrt(1-x^2).
x2=x.*x;


% Calculate efficiently P_l^m (x)/sqrt(1-x^2)^(m/2) - that is, only the
% polynomial part. At least one coefficient is guaranteed to exist - there
% is no null Legendre polynomial.
p=cl*ones(size(x));

for j=l-1:-1:0,
    % Check the exponent of x for current coefficient, px. If it is 0 or 1,
    % just exit the loop
    if(px<2), break; end;
    % If current exponent is >=2, there is a "next" coefficient; multiply p
    % by x2 and add it. Calculate the current coefficient
    cl=-(j+j+2-l-m)*(j+j+1-l-m)/(2*(j+j+1)*(l-j))*cl;
    
    if(maxcf<abs(cl)), maxcf=abs(cl); end;
    %fprintf('Coef is %.16f\n',cl);
    % ...and add to the polynomial
    p=p.*x2 + cl;
    % Decrease the exponent of x - this is the exponent of x corresponding
    % to the newly added coefficient
    px=px-2;
end;
% Estimate the error
err=maxcf*eps;
% fprintf('Coef is %.16f, err %.16f\n',maxcf, err);

% Now we're done adding coefficients. However, if the exponent of x
% corresponding to the last added coefficient is 1 (polynomial is odd),
% multiply the polynomial by x 
if(px==1), p=p.*x; end;

% All that's left is to multiply the whole thing with sqrt(1-x^2)^(m/2). No
% further calculations are needed if m=0.
if(m==0), return; end;

x2=1-x2;
%First, multiply by the integer part of m/2
for j=1:floor(m/2), p=p.*x2; end;
%If m is odd, there is an additional factor sqrt(1-x^2)
if(m~=2*floor(m/2)), p=p.*sqrt(x2); end;

% Finally, the polynomials are not defined for |x|>1. If you do not need
% this behaviour, comment the following line
p(abs(x)>1)=NaN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function to find the coherence of spherical harmonics matrix
% given proposed sampling patterns
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s%%%%%%%%%%%%%%%%
    
function y = simple_objective(x,ang,lm)

ang_proposed=[ang(:,1) x];
[Y,As_Leg]=SH_matrix(ang_proposed,lm,'complex');
y=Tes_Coherence(Y);
end
%% Coherence of a matrix
function val=Tes_Coherence(Basis)
    % absA=abs(A'*A);
    A=Basis;
    normA=zeros(size(A,1),size(A,2));
    for ii=1:size(A,2);
        normA(:,ii)=A(:,ii)/norm(A(:,ii));
    end
    %% Coherence
    InProd=abs(normA'*normA);
    Coher=InProd-eye(size(InProd,1));
    val=max(max(Coher));
end   