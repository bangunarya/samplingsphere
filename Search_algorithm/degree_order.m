function [lmn,lm]=degree_order(L)
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
%% lm for Spherical Harmonics 
l=0:L-1;
 
N=L^2;
lm=zeros(N,2);
idx_beg=1;
 
for o=1:length(l);
    m=-l(o):l(o);
    vec=(combvec(l(o),m)');
    idx=size(vec,1);
    idx_end=idx_beg+idx-1;
    lm(idx_beg:idx_end,:)=vec;
    idx_beg=idx_beg+idx;
     
end
 
 
end