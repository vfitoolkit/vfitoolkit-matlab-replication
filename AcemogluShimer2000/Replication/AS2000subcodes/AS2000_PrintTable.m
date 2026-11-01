function AS2000_PrintTable(TableName,Outs,Paper,fid)
% Prints one of Tables 2-6: replication, with the paper's number in brackets.
% Outs: cell of 4 Out structs (Benchmark, High UI, Long UI, No UI)
% Paper: 18-by-4 matrix of the paper's numbers (NaN where the paper has a dash), or 19-by-4 for Table 6, which has the
% extra row 'Search hours in bad jobs' after 'Work hours in bad jobs'

ColHeaders={'Benchmark','High UI','Long UI','No UI'};
RowNames={'Unemployment benefit b','Benefit duration 1/phi','Tax rate tau','Mean assets',...
    'Fraction of good jobs (%)','Unemployment rate (%)','Unemployment duration','Std. dev. log weekly income',...
    'Work hours in good jobs','Work hours in bad jobs','Search hours by insured','Search hours by uninsured',...
    'Consump. drop following job loss (%)','Std. dev. log annual income',...
    'Average wage','Output','Leisure','Welfare: consumption metric','Search hours in bad jobs'};
% Order the rows are printed in (row 19 only for Table 6)
if size(Paper,1)==19
    RowOrder=[1:10,19,11:18];
else
    RowOrder=1:18;
end

fprintf(fid,'%s \n',TableName);
fprintf(fid,'%-38s',''); for cc=1:4, fprintf(fid,'%20s',ColHeaders{cc}); end; fprintf(fid,'\n');
for ii=1:length(RowOrder)
    rr=RowOrder(ii);
    fprintf(fid,'%-38s',RowNames{rr});
    for cc=1:4
        O=Outs{cc};
        isNoUI=(O.uib==0);
        switch rr
            case 1, x=O.uib;
            case 2, x=1/O.phi;
            case 3, x=O.tau;
            case 4, x=O.meanassets;
            case 5, x=100*O.fracgood;
            case 6, x=100*O.urate;
            case 7, x=O.uduration;
            case 8, x=O.sdlogweeklyincome;
            case 9, x=O.hoursgood;
            case 10, x=O.hoursbad;
            case 11, x=O.hoursinsured;
            case 12
                if isNoUI
                    x=O.hoursallunemployed; % with no UI the insured state is the same as uninsured
                else
                    x=O.hoursuninsured;
                end
            case 13, x=100*O.consdrop;
            case 14, x=O.sdlogannualincome;
            case 15, x=O.avgwage;
            case 16, x=O.output;
            case 17, x=O.leisure;
            case 18, x=O.welfare;
            case 19, x=O.searchbad;
        end
        if isNoUI && any(rr==[2,11])
            x=NaN; % not defined with no UI
        end
        if isNoUI && rr==9 && O.fracgood<1e-6
            x=NaN; % no good jobs
        end
        if isnan(Paper(ii,cc))
            fprintf(fid,'%12.3f (  -   )',x);
        else
            fprintf(fid,'%12.3f (%6.3f)',x,Paper(ii,cc));
        end
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');

end
