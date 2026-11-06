function x=DFJ2010_HRSextract_col(T,name)
% Column 'name' of table T as a double column, or a NaN column if T has no such variable
% (used by DFJ2010_HRSextract for RAND HRS variables that are not asked in every wave).
if ismember(name,T.Properties.VariableNames)
    x=T.(name);
else
    x=nan(height(T),1);
end
end
