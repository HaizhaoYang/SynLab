function peakloc = locatePeaks(txt,prom,blocksize,maxind,method,is_display)

    fft_txt = fftshift(fft2(txt));
    
    h = fspecial('gaussian',15,1);
    fft_txt = log(abs(fft_txt));
%     fft_textb = zeros(size(fft_txt,1)+30,size(fft_txt,2)+30);
%     fft_textb(16:end-15,16:end-15) = fft_txt;
    fft_txt = imfilter(fft_txt,h);
%     fft_txt = fft_txt(16:end-15,16:end-15);
    
    %% find peaks by projections
    if method ==1
        fft_txt1 = max(fft_txt(:,round(end/2):end-20));
        fft_txt2 = max(fft_txt(round(end/2):end-20,:)');
        PeakSig1 = fft_txt1-smooth(fft_txt1,60)';
        PeakSig2 = fft_txt2-smooth(fft_txt2,60)';
    else
        fft_txt1 = sum(fft_txt(:,round(end/2):end-20));
        fft_txt2 = sum(fft_txt(round(end/2):end-20,:),2)';
        PeakSig1 = (smooth(fft_txt1,15)-smooth(fft_txt1,30))/size(fft_txt,1);
        PeakSig2 = (smooth(fft_txt2,15)-smooth(fft_txt2,30))/size(fft_txt,2);
        
    end
    

    % use matlab integrated findpeaks function
    
%     PeakSig1 = smooth(fft_txt1,5)-smooth(fft_txt1,60);
%     PeakSig1 = smooth(fft_txt1,5)-smooth(fft_txt1,30);
    x = 1:length(PeakSig1);
    [pksv,locsv]=findpeaks(PeakSig1,x,'SortStr','descend','MinPeakHeight',0,'MinPeakProminence',prom,'MinPeakDistance',10,'Annotate','extents');
    
    if is_display
        figure, findpeaks(PeakSig1,x,'SortStr','descend','MinPeakHeight',0,'MinPeakProminence',prom,'MinPeakDistance',10,'Annotate','extents');
    end
    
    [pksv,ind] = sort(pksv,'descend');
    locsv = locsv(ind);
    maxind = min(length(locsv),maxind);
    locsv = locsv(1:maxind);

%     PeakSig2 = smooth(fft_txt2,10)-smooth(fft_txt2,60);
%     PeakSig2 = smooth(fft_txt2,5)-smooth(fft_txt2,30);
%     PeakSig2 = fft_txt2-smooth(fft_txt2,60)';
    x = 1:length(PeakSig2);
    [pksh,locsh]=findpeaks(PeakSig2,x,'SortStr','descend','MinPeakHeight',0,'MinPeakProminence',prom,'MinPeakDistance',15,'Annotate','extents');
    
    if is_display
        figure, findpeaks(PeakSig2,x,'SortStr','descend','MinPeakHeight',0,'MinPeakProminence',prom,'MinPeakDistance',15,'Annotate','extents');
    end
    
    [pksh,ind] = sort(pksh,'descend');
    locsh = locsh(ind);
    maxind = min(length(locsh),maxind);
    locsh = locsh(1:maxind);

    % find intersections

    nlines_h = length(locsh);
    nlines_v = length(locsv);
    lines_h = zeros(2*nlines_h,4); 
    lines_v = zeros(2*nlines_v,4);

    for i = 1:nlines_h
        lines_h(2*i-1,:) = [1,floor(blocksize(1)/2)+1+locsh(i),blocksize(2),floor(blocksize(1)/2)+1+locsh(i)];
        lines_h(2*i,:) = [1,floor(blocksize(1)/2)+1-locsh(i),blocksize(2),floor(blocksize(1)/2)+1-locsh(i)];
    end

    for i = 1:nlines_v
        lines_v(2*i-1,:) = [floor(blocksize(2)/2)+1+locsv(i),1,floor(blocksize(2)/2)+1+locsv(i),blocksize(1)];
        lines_v(2*i,:) = [floor(blocksize(2)/2)+1-locsv(i),1,floor(blocksize(2)/2)+1-locsv(i),blocksize(1)];
    end

    lines_v = [lines_v; floor(blocksize(2)/2)+1,1,floor(blocksize(2)/2)+1,blocksize(1)];
    lines_h = [lines_h; 1,floor(blocksize(1)/2)+1,blocksize(2),floor(blocksize(1)/2)+1];

    % calculate intersections
    peakloc = [];
    counter = 0;
    for j = 1:size(lines_h,1)
        for i = 1:size(lines_v,1)

            line_h = lines_h(j,:);
            line_v = lines_v(i,:);

            intersect = InterX([line_h(1),line_h(3);line_h(2),line_h(4)],[line_v(1),line_v(3);line_v(2),line_v(4)]);
            if ~(intersect(2)==floor(blocksize(1)/2)+1 && intersect(1)==floor(blocksize(2)/2)+1)
    %             if ~(abs(intersect(2)-floor(blocksize(2)/2)+1) < lowFq && abs(intersect(1)-floor(blocksize(1)/2)+1) <lowFq)
                counter= counter+1;
                peakloc(counter,:) = [intersect(2),intersect(1)];
            end
        end
    end
end