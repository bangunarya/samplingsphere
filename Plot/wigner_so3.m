%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Function to generate Wigner-D function matrix
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s    
    
function [Wigner_SO3,small_d]=wigner_so3(ang,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Wigner Basis expansion from "Sparse recovery in Wigner-D basis
% expansion 
%       
% Created by Arya Bangun at TI RWTH Aachen
% Last modification: 28.08.2018 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[lmn]= permut_lmn(B);
    function val=eta(m,mu)
        if mu >= m
            val=1;
        else
            val=(-1)^(mu-m);
        end
    end
dmm=zeros(size(ang,1),size(lmn,1));
basis=zeros(size(ang,1),size(lmn,1));
for ii=1:size(lmn,1)
    order=lmn(ii,:);
    miu_plus=abs(order(2)-order(3));
    vau_plus=abs(order(2)+order(3));
    s_plus=order(1)-(miu_plus+vau_plus)/2;
    dmm(:,ii)=sqrt((2*order(1)+1)/(8*pi^2))*eta(order(2),order(3))*sqrt((factorial(s_plus)*factorial(s_plus+miu_plus+vau_plus))/(factorial(s_plus+miu_plus)*factorial(s_plus+vau_plus)))...
        *(sin(ang(:,1)/2).^miu_plus).*(cos(ang(:,1)/2).^vau_plus).*japoly(s_plus,miu_plus,vau_plus,cos(ang(:,1)));
    basis(:,ii)=exp(-1i*order(3)*ang(:,3)).*dmm(:,ii).*exp(-1i*order(2)*ang(:,2)); 
end


small_d=dmm;
Wigner_SO3=basis;
end

function [varargout]=japoly(n,alp,bet,x)
%y=japoly(n,alp,bet,x) computes the Jacobi polynomial of degree n with parameter (alp,bet) at a
%   vector-valued x
% [dy,y]=japoly(n,alp,bet, x) also returns the values of the 1st-order 
%   derivative stored in dy
% See Page 74 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
% 
% Last modified on September 2, 2011    



  apb=alp+bet;     
   
if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=0.5*(alp-bet+(apb+2)*x); return; end;

     polylst=ones(size(x));	
     poly=0.5*(alp-bet+(apb+2)*x);
     
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	  polyn=((a2+a3*x).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
      polylst=poly; poly=polyn;	
   end;
      varargout{1}=polyn; return;
end;

   
if n==0, varargout{2}=ones(size(x)); varargout{1}=zeros(size(x)); return; end;
if n==1, varargout{2}=0.5*(alp-bet+(apb+2)*x); varargout{1}=0.5*(apb+2)*ones(size(x)); return; end;

   polylst=ones(size(x));         pderlst=zeros(size(x));
   poly=0.5*(alp-bet+(apb+2.)*x); pder=0.5*(apb+2.)*ones(size(x));
   
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	  polyn=((a2+a3*x).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
	  pdern=((a2+a3*x).*pder-a4*pderlst+a3*poly)/a1;	  	  
	  polylst=poly; poly=polyn;
	  pderlst=pder; pder=pdern;
   end;
      varargout{2}=polyn;
      varargout{1}=pdern;
   return;

end

function [lmn]= permut_lmn(L)
%% lmn for Wigner Basis
l=0:L-1;
N=L*(2*L-1)*(2*L+1)/3;
lmn=zeros(N,3);
idx_beg=1;
for o=1:length(l);
    m=-l(o):l(o);
    vec=fliplr(combvec(m,m)');
    idx=size(vec,1);
    idx_end=idx_beg+idx-1;
    lmn(idx_beg:idx_end,:)=[l(o)*ones(idx,1) vec];
    idx_beg=idx_beg+idx;
end
end