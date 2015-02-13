clear all;
close all;

%% set up data
load data/synthetic.txt;
2^floor(log2(length(synthetic)))
N = min(2^13,2^floor(log2(length(synthetic))));

x = [0:N-1]/N;
fff = synthetic(1:N);
fff = fff(:).';
f2 = fff/max(fff);

is_real = 1;
alpha = 1; disper = 0.9; delta = 1; beta = 0; NDim = 1; 
ns = RFI_MakeDataAlphaStable(alpha, disper, delta, N, beta, NDim);
ns = ns/max(ns)*15;

if is_real
    ns = real(ns);
    f2 = real(f2);
end

pic = figure;plot(x,f2+ns);axis tight;axis square;xlabel('x');title('Noisy data');

%% set up parameters
res = 4;
NG = N/4;
is_unif = 1;
typeNUFFT = 1;
is_cos = 1;
epsl = 1e-2;
xo = x;
dd = 1.6;

redglobal = 16;
highglobal = N/4;
highband = N/4;
lowband = 100;
for band = 1
    for scnt = 1:3
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% clean signal
        if (1)
            %synchrosqueezed wave packet transform
            red = 16;
            R_high = highglobal;
            R_low = 0;
            if is_real
                freq_range = [R_low R_high];
            else
                freq_range = [-R_high R_high];
            end
            t_sc = 1/2 + scnt/8;
            fff = f2;
            rad = dd + 0.1 -scnt*0;
            [T_f coef kk] = ss_wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,epsl,res,0);
            pic = figure;imagesc([0 1],freq_range,real(T_f));colorbar;title('sswpt clean');axis square;axis xy;
            xlabel('x'); ylabel('k');
            %colormap (1-gray);
            caxis([0,10]);
        end
        
        if (1)
            rad = dd -scnt*0.15;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% signal + Gaussian white noise
            fff = f2 + ns;
            red = 1;
            t_sc = 1/2 + scnt/8;
            if band
                R_low = lowband; R_high = highband;
            else
                R_low = 0; R_high = highglobal;
            end
            ccc = wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
            aaa = wp1_ext(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
            if is_real
                fqscale = [0 highglobal];
            else
                fqscale = [-highglobal highglobal];
            end
            num_grid = ceil((fqscale(2)-fqscale(1))/res);
            h = (fqscale(2)-fqscale(1))/(num_grid-1);
            
            EXT = 10^10;
            ncl = cell(1,red);
            nclp = cell(1,red);
            NG = length(ccc{1}{1}{1});
            EXT = 10^10;
            kk = cell(1,red);
            T_f = zeros(num_grid,NG);
            coef = cell(1,red);
            for cntred = 1:red
                ncl{cntred} = numel(ccc{cntred});
                nclp{cntred} = 0;
                for cnt = 1:ncl{cntred}
                    nclp{cntred} = nclp{cntred} + numel(ccc{cntred}{cnt});
                end
                
                if is_real
                    cc = zeros(nclp{cntred}/2,NG);
                    aa = zeros(nclp{cntred}/2,NG);
                    cnt_nclp = 1;
                    for cnt = 1:ncl{cntred}
                        cc(cnt_nclp,:) = ccc{cntred}{cnt}{1};
                        aa(cnt_nclp,:) = aaa{cntred}{cnt}{1};
                        cnt_nclp = cnt_nclp + 1;
                    end
                else
                    cc = zeros(nclp{cntred}-1,NG);
                    aa = zeros(nclp{cntred}-1,NG);
                    cnt_nclp1 = 1;
                    cnt_nclp2 = 1;
                    for cnt = 1:ncl{cntred}
                        for cnt2 = 1:2
                            if cnt2 == 1
                                cc(cnt_nclp1+nclp{cntred}/2-1,:) = ccc{cntred}{cnt}{cnt2};
                                aa(cnt_nclp1+nclp{cntred}/2-1,:) = aaa{cntred}{cnt}{cnt2};
                                cnt_nclp1 = cnt_nclp1 + 1;
                            else if cnt ~= 1
                                    cc(-cnt_nclp2+nclp{cntred}/2,:) = ccc{cntred}{cnt}{cnt2};
                                    aa(-cnt_nclp2+nclp{cntred}/2,:) = aaa{cntred}{cnt}{cnt2};
                                    cnt_nclp2 = cnt_nclp2 + 1;
                                end
                            end
                        end
                    end
                end
                
                szc = size(cc);
                gud = find(abs(cc)>epsl);
                aa=real(aa(gud)./(2*pi*i*cc(gud)));
                if is_real
                    good = find(aa<=R_high & aa>=R_low);
                else
                    good = find(aa<=R_high & aa>=-R_high);
                end
                
                kk{cntred} = repmat(EXT,szc);
                kk{cntred}(gud(good)) = aa(good);
                
                for b = 1:NG
                    tk = kk{cntred}(:,b);
                    gud = find(tk<EXT);
                    tk = tk(gud);
                    pos = round((tk-fqscale(1))/h)+1;
                    for cnt = 1:length(pos)
                        T_f(pos(cnt),b) = T_f(pos(cnt),b) + abs(cc(gud(cnt),b)).^2;
                    end
                end
            end
            head = sprintf('sswpt energy distribution,red = %d',red);
            %%close all;
            pic = figure;imagesc([0 1],fqscale,real(T_f));colorbar;title(head);axis square; axis xy;%colormap (1-gray);
            xlabel('x'); ylabel('k');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% signal + Gaussian white noise
            fff = f2 + ns;
            red = redglobal;
            t_sc = 1/2 + scnt/8;
            if band
                R_low = lowband; R_high = highband;
            else
                R_low = 0; R_high = highglobal;
            end
            ccc = wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
            aaa = wp1_ext(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
            if is_real
                fqscale = [0 highglobal];
            else
                fqscale = [-highglobal highglobal];
            end
            num_grid = ceil((fqscale(2)-fqscale(1))/res);
            h = (fqscale(2)-fqscale(1))/(num_grid-1);
            
            EXT = 10^10;
            ncl = cell(1,red);
            nclp = cell(1,red);
            NG = length(ccc{1}{1}{1});
            EXT = 10^10;
            kk = cell(1,red);
            T_f = zeros(num_grid,NG);
            coef = cell(1,red);
            for cntred = 1:red
                ncl{cntred} = numel(ccc{cntred});
                nclp{cntred} = 0;
                for cnt = 1:ncl{cntred}
                    nclp{cntred} = nclp{cntred} + numel(ccc{cntred}{cnt});
                end
                
                if is_real
                    cc = zeros(nclp{cntred}/2,NG);
                    aa = zeros(nclp{cntred}/2,NG);
                    cnt_nclp = 1;
                    for cnt = 1:ncl{cntred}
                        cc(cnt_nclp,:) = ccc{cntred}{cnt}{1};
                        aa(cnt_nclp,:) = aaa{cntred}{cnt}{1};
                        cnt_nclp = cnt_nclp + 1;
                    end
                else
                    cc = zeros(nclp{cntred}-1,NG);
                    aa = zeros(nclp{cntred}-1,NG);
                    cnt_nclp1 = 1;
                    cnt_nclp2 = 1;
                    for cnt = 1:ncl{cntred}
                        for cnt2 = 1:2
                            if cnt2 == 1
                                cc(cnt_nclp1+nclp{cntred}/2-1,:) = ccc{cntred}{cnt}{cnt2};
                                aa(cnt_nclp1+nclp{cntred}/2-1,:) = aaa{cntred}{cnt}{cnt2};
                                cnt_nclp1 = cnt_nclp1 + 1;
                            else if cnt ~= 1
                                    cc(-cnt_nclp2+nclp{cntred}/2,:) = ccc{cntred}{cnt}{cnt2};
                                    aa(-cnt_nclp2+nclp{cntred}/2,:) = aaa{cntred}{cnt}{cnt2};
                                    cnt_nclp2 = cnt_nclp2 + 1;
                                end
                            end
                        end
                    end
                end
                
                szc = size(cc);
                gud = find(abs(cc)>epsl);
                aa=real(aa(gud)./(2*pi*i*cc(gud)));
                if is_real
                    good = find(aa<=R_high & aa>=R_low);
                else
                    good = find(aa<=R_high & aa>=-R_high);
                end
                
                kk{cntred} = repmat(EXT,szc);
                kk{cntred}(gud(good)) = aa(good);
                
                for b = 1:NG
                    tk = kk{cntred}(:,b);
                    gud = find(tk<EXT);
                    tk = tk(gud);
                    pos = round((tk-fqscale(1))/h)+1;
                    for cnt = 1:length(pos)
                        T_f(pos(cnt),b) = T_f(pos(cnt),b) + abs(cc(gud(cnt),b)).^2;
                    end
                end
            end
            head = sprintf('sswpt energy distribution,red = %d',red);
            %%close all;
            pic = figure;imagesc([0 1],fqscale,real(T_f));colorbar;title(head);axis square; axis xy;%colormap (1-gray);
            xlabel('x'); ylabel('k');
            % colormap (1-gray);
            caxis([0,1.2+0.3*scnt]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% signal + Gaussian white noise
            if (0)
                fff = f2 + ns;
                red = redglobal;
                t_sc = 1/2 + scnt/8;
                if band
                    R_low = lowband; R_high = highband;
                else
                    R_low = 0; R_high = highglobal;
                end
                ccc = wp1_fwd(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
                aaa = wp1_ext(fff,is_real,is_unif,typeNUFFT,xo,NG,R_high,R_low,rad,is_cos,t_sc,red,0);
                if is_real
                    fqscale = [0 highglobal];
                else
                    fqscale = [-highglobal highglobal];
                end
                num_grid = ceil((fqscale(2)-fqscale(1))/res);
                h = (fqscale(2)-fqscale(1))/(num_grid-1);
                
                EXT = 10^10;
                ncl = cell(1,red);
                nclp = cell(1,red);
                NG = length(ccc{1}{1}{1});
                EXT = 10^10;
                kk = cell(1,red);
                T_f = zeros(num_grid,NG);
                coef = cell(1,red);
                for cntred = 1:red
                    ncl{cntred} = numel(ccc{cntred});
                    nclp{cntred} = 0;
                    for cnt = 1:ncl{cntred}
                        nclp{cntred} = nclp{cntred} + numel(ccc{cntred}{cnt});
                    end
                    
                    if is_real
                        cc = zeros(nclp{cntred}/2,NG);
                        aa = zeros(nclp{cntred}/2,NG);
                        cnt_nclp = 1;
                        for cnt = 1:ncl{cntred}
                            cc(cnt_nclp,:) = ccc{cntred}{cnt}{1};
                            aa(cnt_nclp,:) = aaa{cntred}{cnt}{1};
                            cnt_nclp = cnt_nclp + 1;
                        end
                    else
                        cc = zeros(nclp{cntred}-1,NG);
                        aa = zeros(nclp{cntred}-1,NG);
                        cnt_nclp1 = 1;
                        cnt_nclp2 = 1;
                        for cnt = 1:ncl{cntred}
                            for cnt2 = 1:2
                                if cnt2 == 1
                                    cc(cnt_nclp1+nclp{cntred}/2-1,:) = ccc{cntred}{cnt}{cnt2};
                                    aa(cnt_nclp1+nclp{cntred}/2-1,:) = aaa{cntred}{cnt}{cnt2};
                                    cnt_nclp1 = cnt_nclp1 + 1;
                                else if cnt ~= 1
                                        cc(-cnt_nclp2+nclp{cntred}/2,:) = ccc{cntred}{cnt}{cnt2};
                                        aa(-cnt_nclp2+nclp{cntred}/2,:) = aaa{cntred}{cnt}{cnt2};
                                        cnt_nclp2 = cnt_nclp2 + 1;
                                    end
                                end
                            end
                        end
                    end
                end
                szc = size(cc);
                
                for b = 1:szc(2)
                    gud = find(abs(cc(:,b))>epsl);
                    temp = abs(cc(gud,b));
                    temp = temp(:)';
                    [val ord] = sort(temp,'descend');
                    nb = ceil(length(gud)*0.1);
                    nb = 1;
                    pos = (ord(1:nb));
                    kk = real(sum(aa(gud(pos),b))./sum(cc(gud(pos),b))/(2*pi*i));
                    if is_real
                        %good = find(aa<=R_high & aa>=R_low);
                        if ~(kk<=R_high & kk>=R_low)
                            kk = EXT;
                        end
                    else
                        %good = find(aa<=R_high & aa>=-R_high);
                        if ~(kk<=R_high & kk>=-R_high)
                            kk = EXT;
                        end
                    end
                    if kk<EXT
                        loc = round((kk-fqscale(1))/h)+1;
                        T_f(loc,b) = T_f(loc,b) + sum(abs(cc(gud(pos),b)).^2);
                    end
                end
                head = sprintf('sswpt energy distribution,red = %d,max',red);
                %close all;
                pic = figure;imagesc([0 1],fqscale,real(T_f));colorbar;title(head);axis square; axis xy;%colormap (1-gray);
                xlabel('x'); ylabel('k');
            end
        end
    end
end
