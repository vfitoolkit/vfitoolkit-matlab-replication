function AS2000_LatexTables(SavedOutputDir)
% Writes the LaTeX tables of the write-up to SavedOutputDir/LatexInputs/, from the saved results
% (AS2000_doPart1.mat,...,AS2000_doPart5.mat). Solves nothing.
% Tables 2-6 are replication only (the write-up shows the paper's original table above each one).
% The Section 4.5 table has the paper's numbers in brackets, as the paper only reports them in its text.

if ~exist([SavedOutputDir,'/LatexInputs'],'dir'), mkdir([SavedOutputDir,'/LatexInputs']); end

RowNames={'Unemployment benefit $b$','Benefit duration $1/\phi$','Tax rate $\tau$','Mean assets',...
    'Fraction of good jobs','Unemployment rate','Unemployment duration','Std. dev. of log weekly income',...
    'Work hours in good jobs','Work hours in bad jobs','Search hours by insured','Search hours by uninsured',...
    'Consump. drop following job loss','Std. dev. of log annual income',...
    'Average wage','Output','Leisure','Welfare: consumption metric','Search hours in bad jobs'};
RowFormats={'%.2f','%.1f','%.3f','%.1f','%.1f\\%%','%.1f\\%%','%.1f','%.3f','%.3f','%.3f','%.3f','%.3f',...
    '%.1f\\%%','%.3f','%.3f','%.3f','%.3f','%.3f','%.3f'};
% Blank lines between the blocks of rows, as in the paper's tables
BlockEnds=[4,8,12,14];

for tt=1:5
    load([SavedOutputDir,'/AS2000_doPart',num2str(tt),'.mat'],'Outs')
    % The saved results may hold gpuArrays, bring every field to the cpu
    for cc=1:4
        fnames=fieldnames(Outs{cc});
        for jj=1:length(fnames)
            if isa(Outs{cc}.(fnames{jj}),'gpuArray')
                Outs{cc}.(fnames{jj})=gather(Outs{cc}.(fnames{jj}));
            end
        end
    end
    if tt==5
        RowOrder=[1:10,19,11:18]; % Table 6 has 'Search hours in bad jobs' after 'Work hours in bad jobs'
        BlockEndsT=[4,8,13,15];
    else
        RowOrder=1:18;
        BlockEndsT=BlockEnds;
    end
    fid=fopen([SavedOutputDir,'/LatexInputs/AS2000_Table',num2str(tt+1),'.tex'],'w');
    fprintf(fid,'\\begin{tabular}{lcccc} \n\\hline \n');
    fprintf(fid,'Variable & Column I & Column II & Column III & Column IV \\\\ \n');
    fprintf(fid,' & Benchmark & High UI & Long UI & No UI \\\\ \n\\hline \n');
    for ii=1:length(RowOrder)
        rr=RowOrder(ii);
        fprintf(fid,'%s',RowNames{rr});
        for cc=1:4
            O=Outs{cc};
            isNoUI=(gather(O.uib)==0);
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
            x=double(gather(x));
            if isNoUI && any(rr==[2,11])
                x=NaN; % not defined with no UI
            end
            if isNoUI && rr==9 && gather(O.fracgood)<1e-6
                x=NaN; % no good jobs
            end
            if isnan(x)
                fprintf(fid,' & --');
            else
                fprintf(fid,[' & ',RowFormats{rr}],x);
            end
        end
        fprintf(fid,' \\\\ \n');
        if any(ii==BlockEndsT)
            fprintf(fid,'[2mm] \n');
        end
    end
    fprintf(fid,'\\hline \n\\end{tabular} \n');
    fclose(fid);
end

%% Section 4.5: gains from UI, and welfare including the transition
load([SavedOutputDir,'/AS2000_doPart1.mat'],'Outs')
load([SavedOutputDir,'/AS2000_Section45.mat'],'S45')
% The saved results may hold gpuArrays, bring every field to the cpu
for cc=1:4
    fnames=fieldnames(Outs{cc});
    for jj=1:length(fnames)
        if isa(Outs{cc}.(fnames{jj}),'gpuArray')
            Outs{cc}.(fnames{jj})=gather(Outs{cc}.(fnames{jj}));
        end
    end
end
PaperS45=[17.6, 0.0035, 0.7, 0.6, 2, 0.652; 5.6, 0.0011, 0.5, 0.4, 4, 0.652]; % the paper's text (pg 1215-1216); the ratios are 'twice' and 'about four times'
ColNames={'Benchmark','High UI','Long UI','No UI'};
O1=Outs{1};
UI1=double(gather(O1.UIspending)); Y1=double(gather(O1.output)); W1=double(gather(O1.welfare));
fid=fopen([SavedOutputDir,'/LatexInputs/AS2000_Section45.tex'],'w');
fprintf(fid,'\\begin{tabular}{lcc} \n\\hline \n');
fprintf(fid,' & High UI & Long UI \\\\ \n\\hline \n');
Rows={'Change in UI spending (\\%%)','Change in UI spending per worker','Change in output (\\%%)','Change in welfare (\\%%)',...
    'Change in output / change in UI spending','Welfare including the transition'};
Formats={'%.1f','%.4f','%.2f','%.2f','%.2f','%.4f'};
PaperFormats={'%.1f','%.4f','%.1f','%.1f','%.0f','%.3f'};
for rr=1:6
    fprintf(fid,Rows{rr});
    for cc=2:3
        O=Outs{cc};
        UI=double(gather(O.UIspending)); Y=double(gather(O.output)); W=double(gather(O.welfare));
        switch rr
            case 1, x=100*(UI/UI1-1);
            case 2, x=UI-UI1;
            case 3, x=100*(Y/Y1-1);
            case 4, x=100*(W/W1-1);
            case 5, x=(Y-Y1)/(UI-UI1);
            case 6, x=S45.welfaretransition(cc);
        end
        fprintf(fid,[' & ',Formats{rr},' (',PaperFormats{rr},')'],x,PaperS45(cc-1,rr));
    end
    fprintf(fid,' \\\\ \n');
end
fprintf(fid,'\\hline \n\\end{tabular} \n');
fclose(fid);

% Welfare in the benchmark, for the text
fid=fopen([SavedOutputDir,'/LatexInputs/AS2000_Section45_benchmarkwelfare.tex'],'w');
fprintf(fid,'%.4f',W1);
fclose(fid);

end
