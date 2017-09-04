function [mloc1, mloc2] = sift_match(img1, img2, thresh)
    [des1, loc1] = sift(img1);
    [des2, loc2] = sift(img2);
    N1 = size(des1,1);
    match_table = zeros(1,N1);
    for i = 1:N1
        dotprods = des1(i,:) * des2';
        [sort_dis,idx] = sort(acos(dotprods));
        if(sort_dis(1) < thresh*sort_dis(2))
            match_table(i) = idx(1);
        else
            match_table(i) = 0;
        end
    end
    num = sum(match_table > 0);
    
    idx1 = find(match_table);
    idx2 = match_table(idx1);
    x1 = loc1(idx1,2);
    x2 = loc2(idx2,2);
    y1 = loc1(idx1,1);
    y2 = loc2(idx2,1);
    mloc1 = [x1,y1];
    mloc2 = [x2,y2];
end

