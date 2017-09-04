function imgout = img_stitch(img1, img2)
    %transform img2 to img1, img1 is steady
    [loc1, loc2] = sift_match(img1, img2, 0.6);
    H = ransac(loc2, loc1);
    tform = maketform('projective',H');
    img21 = imtransform(img2,tform);
    [M1,N1,~] = size(img1);
    [M2,N2,~] = size(img2);
    [M3,N3,~] = size(img21);
    pt1 = [1 1 1;1 M2 1;N2 1 1;N2 M2 1];
    pt2 = H*pt1';
    
    pt2 = pt2(1:2,:)' ./ repmat(pt2(3,:)',1,2);
    x = pt2(:,1);
    y = pt2(:,2);

   up = round(min(y));
    Yoffset = 0;

    if up <= 0
        Yoffset = -up+1;
        up = 1;
    end

    left = round(min(x));
    Xoffset = 0;

    if left<=0
        Xoffset = -left+1;
        left = 1;
    end
    
    imgout(up:up+M3-1,left:left+N3-1,:) = img21;
    
    x = 1;
    y = 1;
    for row = Yoffset+1:Yoffset+M1
        y = 1;
        for col = Xoffset+1:Xoffset+N1
            if(sum(img1(x,y,:))~=0)
                imgout(row,col,:) = img1(x,y,:);
            end
            y = y+1;
        end
        x = x+1;
    end
    
    %imgout(Xoffset+1:Xoffset+M1,Yoffset+1:Yoffset+N1,:) = img1;
    imshow(imgout);
end

