
%%  chi square test  for image data
function [VarP,  ChiP, VarC, ChiC] =Var_Chisquare_Gray(P,C)  
    [M,N]=size(P);   
    Ph=M*N/256;
 
    P=uint16(P);  % P=Plain image
    Pim=imhist(P(:));
    VarP=sum((Pim-Ph).^2)/256;
    ChiP=sum((Pim-Ph).^2)/(Ph);
    
    
    C=uint8(C);  % C=cipher image
    Cim=(imhist(C(:)));
    VarC=sum((Cim-Ph).^2)/256;
    ChiC=sum((Cim-Ph).^2)/(Ph);
end 