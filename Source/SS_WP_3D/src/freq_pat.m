function pat = freq_pat(H,rad,R_low,R_high)
if nargin <2, rad = 4; end;
    
    if(H<16)    error('H needs to be at least 16');  end

        if(1)
            tr = rad;
            tl = [1 1 1 1];
            cnt = 4*rad;
            idx = 1;
        end
        
        pat{1} = struct('wid',tr,'lst',tl,'ed',cnt);
        while(cnt<H)
            old = cnt;
            pat{idx}.lst = [pat{idx}.lst 1 1];
            pat{idx}.ed = pat{idx}.wid*2 + pat{idx}.ed;
            cnt = cnt + 2*tr;
            idx = idx+1;
            rad = 2*rad;
            trg = min(4*old,H);
            tl = [zeros(1,cnt/rad), ones(1,(trg-cnt)/rad)];
            tr = rad;
            pat{idx} = struct('wid',tr,'lst',tl,'ed',trg);
            cnt = trg;
        end
        
        %cut the pat
        cnt = 1;
        while(pat{cnt}.ed<=R_low)
            cnt = cnt + 1;
        end
        cnt2 = length(pat);
        while(pat{cnt2}.ed>=R_high)
            cnt2 = cnt2 -1;
        end
        if cnt2+1<length(pat)
            pat = pat(1:cnt2+1);
        end
        if cnt-1>1
            pat = pat(cnt-1:end);
        end
  
  