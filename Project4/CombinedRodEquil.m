function T = CombinedRodEquil(x,coeff)
%COMBINEDRODEQUIL Summary of this function goes here
%   Detailed explanation goes here
if length(coeff)==3
    coeffAluminum=coeff([1,3]);
    coeffCopper=coeff([2,3]);
else 
    coeffAluminum=coeff;
    coeffCopper=coeff;
end
T(:,1)=UninsulatedRodEquil(x,coeffAluminum,'Aluminum');
T(:,2)=UninsulatedRodEquil(x,coeffCopper,'Copper');
end

