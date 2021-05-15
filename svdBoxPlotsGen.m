function [BxpltPvalue,pkvalue] = svdBoxPlotsGen(File,NumCat,Vrows,pvcut,Cat,CatL)
% U, S, V are components of SVD
% Vrows rows in V starting @ 1 that function should run through
% NumCat is number of categories being tested
% Cat is a matrix of form [('category name'),size,...]
% function should display boxplots and have output of pvalues for each row
% of V chosen to run through
data = csvread(File,4,1);
[U,S,V] = svd(data,0);
dS = diag(S);
V = V';
j = 1;
pkvalue = 0;

for k = 1:Vrows
    if NumCat > 2
        ind1 = V(k,1:(CatL(1)-1));
        ind2 = V(k,CatL(1):(CatL(1)+CatL(2)-1));
        ind3 = V(k,(CatL(1)+CatL(2)):(CatL(1)+CatL(2)+CatL(3)-1));
        ind4 = V(k,(CatL(1)+CatL(3)+CatL(2)):(CatL(1)+CatL(2)+CatL(3)+CatL(4)-1));
        iData = [ind1';ind2';ind3';ind4'];

        i1 = repmat(Cat(1,1),length(ind1),1);
        i2 = repmat(Cat(1,2),length(ind2),1);
        i3 = repmat(Cat(1,3),length(ind3),1);
        i4 = repmat(Cat(1,4),length(ind4),1);
        iCat = [i1;i2;i3;i4];
        
        BxpltPvalue(k) = kruskalwallis(iData,iCat,'off');
        
        if BxpltPvalue(k) < pvcut
            figure(k)
            boxplot(iData,iCat)
            pkvalue(1,j) = k;
            pkvalue(2,j) = BxpltPvalue(k);
            pkvalue(3,j) = dS(k);
            j = j+1;
        end
       
    else
        ib1 = V(k,1:(CatL(1)-1));
        ib2 = V(k,CatL(1):(CatL(1)+CatL(2)-1));
       
        iData = [ib1';ib2'];

        i1 = repmat(Cat(1,1),length(ib1),1);
        i2 = repmat(Cat(1,2),length(ib2),1);
        iCat = [i1;i2];
        
        BxpltPvalue(k) = ranksum(ib1,ib2);
        
        
        if BxpltPvalue(k) < pvcut
            figure(k)
            boxplot(iData,iCat)
            pkvalue(1,j) = k;
            pkvalue(2,j) = BxpltPvalue(k);
            pkvalue(3,j) = dS(k);
            j = j+1;
        end
       
    end
end
