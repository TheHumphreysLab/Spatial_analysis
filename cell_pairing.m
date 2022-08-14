clear all
close all

week6_stat = readtable('~/Desktop/new_coord_week6.csv');
CellList = table2struct(week6_stat);
CellTypeID_new = [CellList.cell_index];
PairCounts = zeros(max(CellTypeID_new),max(CellTypeID_new));
Sum = 0;
Sum_nonself = 0; 
NicheRadius = 200;
for i = 1:length(CellList)-1
    for j = i+1:length(CellList)
        x1 = CellList(i).x;
        y1 = CellList(i).y;
        x2 = CellList(j).x;
        y2 = CellList(j).y;
        if ((x1-x2)^2 + (y1-y2)^2)^0.5 < NicheRadius
            Type1 = CellList(i).cell_index;
            Type2 = CellList(j).cell_index;
            if Type1 ~= Type2
            	PairCounts(Type1,Type2) = PairCounts(Type1,Type2)+1;
                PairCounts(Type2,Type1) = PairCounts(Type2,Type1)+1;
                Sum_nonself = Sum_nonself+1;
            else
                PairCounts(Type1,Type2) = PairCounts(Type1,Type2)+1;
            end
            Sum = Sum+1;
        end            
    end
end
NicheEnrichment_nonself = PairCounts/Sum_nonself;
NicheEnrichment = PairCounts/Sum;

PairCounts0 = zeros(max(CellTypeID_new),max(CellTypeID_new));
Sum0 = 0;
Sum0_nonself = 0; 
for i = 1:length(CellList)-1
    for j = i+1:length(CellList)
        x1 = CellList(i).x;
        y1 = CellList(i).y;
        x2 = CellList(j).x;
        y2 = CellList(j).y;
        Type1 = CellList(i).cell_index;
        Type2 = CellList(j).cell_index;
        if Type1 ~= Type2
            PairCounts0(Type1,Type2) = PairCounts0(Type1,Type2)+1;
            PairCounts0(Type2,Type1) = PairCounts0(Type2,Type1)+1;
            Sum0_nonself = Sum0_nonself+1;
        else
            PairCounts0(Type1,Type2) = PairCounts0(Type1,Type2)+1;
        end
        Sum0 = Sum0+1;
    end
end
NicheEnrichment0_nonself = PairCounts0/Sum0_nonself;
NicheEnrichment0 = PairCounts0/Sum0;
NicheEnrichment = NicheEnrichment./NicheEnrichment0;
NicheEnrichment_nonself = NicheEnrichment_nonself./NicheEnrichment0_nonself;
for i = 1:max(CellTypeID_new)
    NicheEnrichment_nonself(i,i) = 1;
end

