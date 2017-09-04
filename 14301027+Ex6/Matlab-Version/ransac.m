function H = ransac(pt1, pt2)
    N = size(pt1,1);
    max_inline = 0;
    for iter = 1:200
        random_index = randperm(N);
        ori = pt1(random_index(1:round(N/5)),:);
        tar = pt2(random_index(1:round(N/5)),:);
        trans = get_eigen(ori,tar);
        temp = trans * [pt1,ones(N,1)]';
        temp = temp(1:2,:)'./repmat(temp(3,:)',1,2);
        dis = sum((temp-pt2).*(temp-pt2),2);
        inline = sum(dis < 4);
        if inline > max_inline
            max_inline = inline;
            H = trans;
        end
    end
end


