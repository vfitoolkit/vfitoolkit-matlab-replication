function [T,varnames,info]=DFJ2010_ReadStata114(filename,wantvars,verbose)
% Read selected variables from a Stata format-114/115 (Stata 10-12) .dta file.
% The RAND HRS Longitudinal, Detailed Imputations, Exit and Finder files are all
% format 114, which is a fixed-width binary record layout, so we can read only the
% columns we need without loading the (1.7GB) file into memory.
%
% Inputs:
%   filename  path to .dta
%   wantvars  cellstr of variable names to read (case-insensitive). Any not present
%             in the file are reported in info.notfound and skipped.
%             Use {} to read only the header (T is empty, varnames is the full list).
%   verbose   (optional, default 1) print progress
% Outputs:
%   T         table, one column per found variable (numeric double; Stata missing
%             codes ., .a, ..., .z become NaN; string variables become string arrays)
%   varnames  cellstr of ALL variable names in the file (in file order)
%   info      struct: nvar, nobs, reclen, byteorder, datalabel, notfound, varlabels
%
% Stata 114 layout (see Stata's -help dta- for release 10/11):
%   header:  ds_format(1) byteorder(1) filetype(1) unused(1) nvar(int16) nobs(int32)
%            data_label(81) time_stamp(18)
%   descriptors: typlist(nvar bytes) varlist(33*nvar) srtlist(2*(nvar+1)) fmtlist(49*nvar)
%            lbllist(33*nvar)
%   variable labels: 81*nvar
%   expansion fields: repeated {dtype(1) len(int32) contents(len)} until dtype=0,len=0
%   data: nobs records, each the concatenation of the variables in typlist order
% typlist codes: 251 byte(1), 252 int(2), 253 long(4), 254 float(4), 255 double(8),
%            1..244 str of that length.

if nargin<3
    verbose=1;
end

fid=fopen(filename,'r');
if fid<0
    error('DFJ2010_ReadStata114: cannot open %s',filename);
end
cleanupObj=onCleanup(@() fclose(fid));

ds_format=fread(fid,1,'uint8');
if ~(ds_format==114 || ds_format==115)
    error('DFJ2010_ReadStata114: file is Stata format %d, this reader handles 114/115 only',ds_format);
end
byteorder=fread(fid,1,'uint8'); % 1=HILO (big endian), 2=LOHI (little endian)
if byteorder==2
    mf='ieee-le';
else
    mf='ieee-be';
end
fread(fid,1,'uint8'); % filetype (=1)
fread(fid,1,'uint8'); % unused
nvar=fread(fid,1,'int16',0,mf);
nobs=fread(fid,1,'int32',0,mf);
datalabel=deblank(char(fread(fid,81,'uint8')'));
fread(fid,18,'uint8'); % time stamp

typlist=fread(fid,nvar,'uint8');
varnames=cell(nvar,1);
raw=fread(fid,[33,nvar],'uint8');
for ii=1:nvar
    s=char(raw(:,ii)');
    varnames{ii}=s(1:find(s==0,1)-1);
end
fread(fid,2*(nvar+1),'uint8'); % srtlist
fread(fid,49*nvar,'uint8'); % fmtlist
fread(fid,33*nvar,'uint8'); % lbllist
varlabels=cell(nvar,1);
raw=fread(fid,[81,nvar],'uint8');
for ii=1:nvar
    s=char(raw(:,ii)');
    z=find(s==0,1);
    if isempty(z)
        varlabels{ii}=s;
    else
        varlabels{ii}=s(1:z-1);
    end
end
% expansion fields
while true
    dtype=fread(fid,1,'uint8');
    len=fread(fid,1,'int32',0,mf);
    if isempty(dtype) || (dtype==0 && len==0)
        break
    end
    fseek(fid,len,'cof');
end
datastart=ftell(fid);

% byte width of each variable and its offset within a record
width=zeros(nvar,1);
width(typlist==251)=1;
width(typlist==252)=2;
width(typlist==253)=4;
width(typlist==254)=4;
width(typlist==255)=8;
isstr=(typlist<=244);
width(isstr)=typlist(isstr);
offset=[0;cumsum(width(1:end-1))];
reclen=sum(width);

info.nvar=nvar;
info.nobs=nobs;
info.reclen=reclen;
info.byteorder=byteorder;
info.datalabel=datalabel;
info.varlabels=varlabels;
info.notfound={};

if verbose
    fprintf('DFJ2010_ReadStata114: %s\n   %d variables, %d observations, record length %d bytes\n',filename,nvar,nobs,reclen);
end

T=table();
if isempty(wantvars)
    return
end

% resolve the requested variables
[found,idx]=ismember(lower(wantvars(:)),lower(varnames));
info.notfound=wantvars(~found);
if any(~found) && verbose
    fprintf('   NOT in file (skipped): %s\n',strjoin(wantvars(~found),' '));
end
idx=idx(found);
names=wantvars(found); % columns are named as requested (the RAND files store names in lowercase)
names=names(:)';
nsel=numel(idx);

% read the data in blocks of records, keep only the selected columns
blocksize=max(1,floor(200e6/reclen)); % ~200MB per block
cols=cell(nsel,1);
for kk=1:nsel
    if isstr(idx(kk))
        cols{kk}=strings(nobs,1);
    else
        cols{kk}=zeros(nobs,1);
    end
end
fseek(fid,datastart,'bof');
done=0;
while done<nobs
    nb=min(blocksize,nobs-done);
    block=fread(fid,[reclen,nb],'*uint8');
    if size(block,2)<nb
        error('DFJ2010_ReadStata114: file ended after %d of %d records',done+size(block,2),nobs);
    end
    for kk=1:nsel
        v=idx(kk);
        b=block(offset(v)+1:offset(v)+width(v),:); % width x nb
        t=typlist(v);
        if isstr(v)
            c=cellstr(char(b'));
            cols{kk}(done+1:done+nb)=string(c);
        else
            if byteorder==1 % big endian: reverse bytes
                b=flipud(b);
            end
            switch t
                case 251
                    x=double(typecast(b(:),'int8'));
                    x(x>100)=NaN; % . is 101, .a=102, ..., .z=127
                case 252
                    x=double(typecast(b(:),'int16'));
                    x(x>32740)=NaN;
                case 253
                    x=double(typecast(b(:),'int32'));
                    x(x>2147483620)=NaN;
                case 254
                    x=double(typecast(b(:),'single'));
                    x(x>1.701e38)=NaN;
                case 255
                    x=typecast(b(:),'double');
                    x(x>8.988e307)=NaN;
            end
            cols{kk}(done+1:done+nb)=x;
        end
    end
    done=done+nb;
    if verbose
        fprintf('   read %d of %d records\n',done,nobs);
    end
end

T=table(cols{:},'VariableNames',names);

end
