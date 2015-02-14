function brewermap_read(man)
% Read the ColorBrewer Excel file and print formatted code into "brewermap.m".
%
% (c) 2014 Stephen Cobeldick
%
% Due to license restrictions it is not permitted to offer the ColorBrewer
% schemes on MATLAB exchange (ColorBrewer uses an Apache license). This
% function imports an excel spreadsheet available from the ColorBewer website
% and prints correctly formatted code directly into the Mfile "brewermap.m".
%
% ### Instructions ###
%
% 1. Save both "brewermap.m" and "brewermap_read.m" in the current folder.
% 2. Download and save in the current folder the following excel spreadsheet:
%    http://www.personal.psu.edu/cab38/ColorBrewer/ColorBrewer_all_schemes_RGBonly3.XLS
% 3. Run "brewermap_read.m". It will insert the correct code into "brewermap.m".
% 4. Check that calling >>brewermap('demo') creates a figure of all colorschemes.
% 5. Finished! You are ready to use the full capabilities of "brewermap.m".
%
% If all steps are complete, then then following files may be deleted:
% - brewermap_read.m
% - brewermap_temp.m
% - ColorBrewer_all_schemes_RGBonly3.XLS
%
% ### Browsable Download ###
%
% Alternatively you can browse to the required spreadsheet:
% 1. Open http://colorbrewer.org/
% 2. Click 'downloads'
% 3. Click 'Download Spreadsheet'
% 4. Save "ColorBrewer_all_schemes_RGBonly3.XLS" in the current folder.
%
% ### If All Else Fails ###
%
% You can also print the required code in the command window by calling
% >>brewermap_read(true). Copy the printed code from the command window
% and paste it into "brewermap.m", replacing the two lines starting with
% 'XXXXX' and 'YYYYY', according to the printed information headings.
%
% As an absolute last resort you can write the switch statement to replace
% the line starting 'XXXXX', where each switch case follows this template:
%    case 'Name' % ColorBrewer colorscheme name
%        rgb = [R1,G1,B1;R2,G2,B2;...]; % RGB 0<=val<=255, one row per node.
%        typ = 'Type'; % 'Diverging'/'Qualitative'/'Sequential'
% And the cell array listing each colorscheme name, to replace line 'YYYYY':
%seq = {'Name1';'Name2';...}; % ColorBrewer colorscheme names
%

if nargin>0&&man % Manual (print in Command Window)
    fid = 1;
    nwl = sprintf('\n');
else % Automatic (print in "brewermap.m")
    nwl = sprintf('\r\n');
    nam = 'brewermap';
    % Read code from "brewermap.m":
    fid = fopen([nam,'.m']);
    str = textscan(fid, '%s', 'Delimiter','\n', 'WhiteSpace','');
    fclose(fid);
    str = str{1};
    % Check if replacement markers exists in code:
    idx = strncmp('XXXXX',str,5);
    idy = strncmp('YYYYY',str,5);
    switch sum(idx)
        case 0
            fprintf('The Mfile <%s.m> is already complete.\n',nam)
            return
        case 1
            % continue to replace line...
        otherwise
            error('Mfile <%s.m> has too many lines starting "XXXXX".',nam)
    end
end
%
% Import Excel spreadsheet into MATLAB:
s = warning('off','MATLAB:xlsread:Mode');
[num,txt] = xlsread('ColorBrewer_all_schemes_RGBonly3.XLS','','','basic');
warning(s)
%
% Separate header from data:
hdr = @(s) find(strcmpi(txt(1,:),s));
txt(1,:) = [];
num(1,:) = [];
% Remove trailing rows:
txt(cellfun('isempty',txt(:,hdr('ColorLetter'))),:) = [];
% Locate indices belonging to one scheme:
ids = find([~cellfun('isempty',txt(:,hdr('SchemeType')));true]);
% Sort by scheme type, then name:
[srt,idr] = sortrows(txt(ids(1:end-1),[hdr('SchemeType'),hdr('ColorName')]));
% Make list of all colorschemes:
seq = sprintf('''%s'';',srt{:,2});
seq = ['seq = {',seq(1:end-1),'};'];
%
if fid>=3
    % Create new token list
    switch sum(idy)
        case 0
            % no line to replace
        case 1
            str{idy} = seq;
        otherwise
            error('Mfile <%s.m> has too many lines starting "YYYYY".',nam)
    end
    % Create temporary file (in case of an error), print leading code:
    fid = fopen([nam,'_temp.m'],'w');
    idz = find(idx);
    for n = 1:idz-1
        fprintf(fid,'%s%s', str{n}, nwl);
    end
else
    fprintf('YYYYY line replacement:\n\n')
    fprintf('%s\n\n',seq)
    fprintf('XXXXX line replacement:\n\n')
end
%
% For each scheme print corectly formatted code:
for n = idr.'
    tmp = txt{ids(n),hdr('ColorName')};
    fprintf(fid,'    case ''%s'' %% %s%s        rgb = [', lower(tmp), tmp, nwl);
    idq = ids(n):ids(n+1)-1;
    [~,idu] = unique(txt(idq,hdr('ColorLetter')));
    idc = [hdr('R'),hdr('G'),hdr('B')];
    fprintf(fid,'%.0f,%.0f,%.0f;', num(idq(idu(1:end-1)),idc).');
    fprintf(fid,'%.0f,%.0f,%.0f];%s        typ = ''%s'';%s',...
        num(idq(idu(end)),idc), nwl, txt{ids(n),hdr('SchemeType')}, nwl);
end
%
if fid>=3
    % Print trailing code:
    for n = idz+1:numel(str)-1
        fprintf(fid,'%s%s', str{n}, nwl);
    end
    fprintf(fid,'%s',str{end});
    fclose(fid);
    % Replace existing "brewermap.m" with new Mfile:
    if movefile([nam,'_temp.m'],[nam,'.m'])
        fprintf('The Mfile <%s.m> is now ready to use!\n',nam)
    else
        fprintf('Could not replace <%s.m>. Ensure Mfile is writable.\n',nam)
    end
end
%
end
%----------------------------------------------------------------------END:brewermap_read