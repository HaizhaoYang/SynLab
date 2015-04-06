function BD = findGB(TTEng_1st,TTEng_2nd,algorithm)
if algorithm == 1
    BD = 1./(TTEng_1st-TTEng_2nd+1);
else
    [N1,N2,num_wave] = size(TTEng_1st);
    temp = zeros(N1,N2);
    for cnt = 1:num_wave
        temp = temp + TTEng_1st(:,:,cnt) + TTEng_2nd(:,:,cnt);
    end
    W = zeros(N1,N2,num_wave);
    for cnt = 1:num_wave
        W(:,:,cnt) = (TTEng_1st(:,:,cnt)+TTEng_2nd(:,:,cnt))./temp;
    end
    BD = zeros(N1,N2);
    for cnt = 1:num_wave
        BD = BD + W(:,:,cnt)./(TTEng_1st(:,:,cnt)-TTEng_2nd(:,:,cnt)+1);
    end
end