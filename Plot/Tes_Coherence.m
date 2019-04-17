%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check mutual coherence of a matrix
%
% Created by Arya Bangun at IHF RWTH Aachen
% Last modification: 28.08.2014 by Arya Bangun
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%s
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
    
 