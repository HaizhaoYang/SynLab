function dat = wp3_setup(Nx,Ny,Nz,SPx,SPy,SPz,rad,R_low,R_high)
% This code generates the partition of time-frequency domain.
%
% Input:
% [Nx Ny Nz] is the size of the image
% [SPx SPy SPz] is the size of the grid in the space
% rad is the smallest size of the supports of wave packets
% [R_low R_high]         The range of interested spectrum
%
% Output:
% dat stores the information of the grid of wave packets
%
% By Haizhao Yang
    N = max([Nx,Ny,Nz]);
    pat = freq_pat(N/2,rad,R_low,R_high);
    
    if(min([SPx,SPy,SPz])<pat{end}.wid * 2)
        error('Wrong, the size of sampling grid in space is too small');
    end
    
    %dat = cell(numel(pat),1);
    cnt = 0;
    dat = cell(100000,1);
    
    %final width, the same sampling density for all bumps
    
    for ell=1:numel(pat)
        B = pat{ell}.wid;    lst = pat{ell}.lst;
        nw = numel(lst);
        nz = sum(lst==0);
        for I=0:nw-1
            for J=0:nw-1
                for K=0:nw-1
                    if(lst(I+1)==0 && lst(J+1)==0 && lst(K+1)==0)
                        %do nothing
                    else
                        Ist = I*B;      Jst = J*B;       Kst = K*B;%starting position in freq
                        Ied = (I+1)*B;      Jed = (J+1)*B;      Ked = (K+1)*B;
                        if(mod(I,2)==0)
                            Ifm = Ist-2/3*B;        Ito = Ist+4/3*B;                      Ict = (Ifm+Ito)/2;
                        else
                            Ifm = Ist-1/3*B;        Ito = Ist+5/3*B;                      Ict = (Ifm+Ito)/2;
                        end
                        if(mod(J,2)==0)
                            Jfm = Jst-2/3*B;        Jto = Jst+4/3*B;                      Jct = (Jfm+Jto)/2;
                        else
                        Jfm = Jst-1/3*B;        Jto = Jst+5/3*B;                      Jct = (Jfm+Jto)/2;
                        end
                        if(mod(K,2)==0)
                            Kfm = Kst-2/3*B;        Kto = Kst+4/3*B;                      Kct = (Kfm+Kto)/2;
                        else
                        Kfm = Kst-1/3*B;        Kto = Kst+5/3*B;                      Kct = (Kfm+Kto)/2;
                        end
                        %all 8 bumps
                        for id=0:1
                            if(id==0)
                                Idx = [ceil(Ifm):floor(Ito)];      Icf = kf_rt(Idx/B*pi, I);                          Ict = Ict;
                            else
                                Idx = [ceil(-Ito):floor(-Ifm)];      Icf = kf_lf(Idx/B*pi, I);                          Ict = -Ict;
                                tmp = Ist;                            Ist = -Ied;                            Ied = -tmp;
                            end
                            for jd=0:1
                                if(jd==0)
                                    Jdx = [ceil(Jfm):floor(Jto)];      Jcf = kf_rt(Jdx/B*pi, J);                              Jct = Jct;
                                else
                                    Jdx = [ceil(-Jto):floor(-Jfm)];      Jcf = kf_lf(Jdx/B*pi, J);                              Jct = -Jct;
                                    tmp = Jst;                                Jst = -Jed;                            Jed = -tmp;
                                end
                                for kd=0:1  
                                    if(kd==0)
                                        Kdx = [ceil(Kfm):floor(Kto)];      Kcf = kf_rt(Kdx/B*pi, K);                              Kct = Kct;
                                    else
                                        Kdx = [ceil(-Kto):floor(-Kfm)];      Kcf = kf_lf(Kdx/B*pi, K);                              Kct = -Kct;
                                        tmp = Kst;                                Kst = -Ked;                            Ked = -tmp;
                                    end
                                    tmp = abs(kron3(Icf,Jcf,Kcf));
                                    %add info
                                    cur = struct('rng',0,'ctr',0,'k1gd',0,'k2gd',0,'k3gd',0,'wgt',0,'cgw',0,'ex1',0,'ex2',0,'ex3',0);
                                    cur.rng = [Ist Jst Kst Ied Jed Ked];
                                    cur.ctr = [Ict Jct Kct];
                                    cur.k1gd = Idx;           
                                    cur.k2gd = Jdx;
                                    cur.k3gd = Kdx;
                                    [t1,t2,t3] = ndgrid(cur.k1gd,cur.k2gd,cur.k3gd);
                                    cur.wgt = tmp;
                                    cur.cgw = tmp;
                                    cur.ex1 = (-2*pi*i*t1).*cur.wgt;
                                    cur.ex2 = (-2*pi*i*t2).*cur.wgt;
                                    cur.ex3 = (-2*pi*i*t3).*cur.wgt;
                                    dat{cnt+1} = cur;
                                    cnt=cnt+1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    dat = dat(1:cnt);
end

