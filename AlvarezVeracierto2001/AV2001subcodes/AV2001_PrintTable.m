function AV2001_PrintTable(TableName,ColHeaders,Outs,RowNames,fid,latex,Base)
% Prints one of the AV2001 tables. Outs is a cell array of Out structs from
% AV2001_StationaryEqm; the index rows are relative to Base (=100), which is
% the first column unless given (Table 1's U.S. column comes before its base). fid is an open file id (or 1 for the screen).
% latex=0: aligned text (the diary). latex=1: a LaTeX tabular for
% ./SavedOutput/LatexInputs/ (just the tabular; the write-up supplies the
% table environment and caption, TableName is not printed).
%
% Rates are in percent; index rows are relative to the first column.
% Welfare: proportionate increase in permanent consumption that makes average
% lifetime utility in the laissez-faire economy equal to that in the economy
% under consideration; under log utility this is exp((W-W_LF)/Lbar)-1.

nCol=length(Outs);
if nargin<7
    Base=Outs{1};
end
if latex==1
    fprintf(fid,'\\begin{tabular}{l%s} \\hline \n',repmat('c',1,nCol));
    fprintf(fid,'Months of wages'); for cc=1:nCol, fprintf(fid,' & %s',ColHeaders{cc}); end; fprintf(fid,' \\\\ \\hline \n');
    pct='%.1f\\%%'; idx='%.1f'; sep=' & '; eol=' \\\\ \n';
else
    fprintf(fid,'%s \n',TableName);
    fprintf(fid,'%-32s',''); for cc=1:nCol, fprintf(fid,'%14s',ColHeaders{cc}); end; fprintf(fid,'\n');
    pct='%13.1f%%'; idx='%14.1f'; sep=''; eol='\n';
end
for rr=1:length(RowNames)
    if latex==1
        fprintf(fid,'%s',RowNames{rr});
    else
        fprintf(fid,'%-32s',RowNames{rr});
    end
    for cc=1:nCol
        O=Outs{cc};
        fprintf(fid,sep);
        switch RowNames{rr}
            case 'Job-finding rate'
                fprintf(fid,pct,100*O.jobfindingrate);
            case 'Layoff rate'
                fprintf(fid,pct,100*O.layoffrate);
            case 'Unemployment rate'
                fprintf(fid,pct,100*O.unemploymentrate);
            case 'Interest rate'
                fprintf(fid,pct,100*O.annualinterestrate);
            case 'Employment'
                fprintf(fid,idx,100*O.N_hh/Base.N_hh);
            case 'Output'
                fprintf(fid,idx,100*O.Ymeasured/Base.Ymeasured); % measured output (excludes investment in new establishments, as in the paper)
            case 'Capital'
                fprintf(fid,idx,100*O.K/Base.K);
            case 'Assets'
                fprintf(fid,idx,100*O.S/Base.S);
            case 'Consumption'
                fprintf(fid,idx,100*O.C/Base.C);
            case 'S.D. consumption'
                fprintf(fid,idx,100*O.sdC/Base.sdC);
            case 'Non-search time'
                fprintf(fid,idx,100*O.nonsearchtime/Base.nonsearchtime);
            case 'Welfare: no search disutility'
                fprintf(fid,idx,100*exp((O.W_noalpha-Base.W_noalpha)/Base.Lbar));
            case 'Welfare'
                fprintf(fid,idx,100*exp((O.W-Base.W)/Base.Lbar));
            otherwise
                error('AV2001_PrintTable: unknown row %s',RowNames{rr})
        end
    end
    fprintf(fid,eol);
end
if latex==1
    fprintf(fid,'\\hline \n\\end{tabular} \n');
else
    fprintf(fid,'\n');
end

end
