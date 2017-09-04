function H = get_affine(pt1, pt2)
    N = size(pt1,1);
    p1 = [pt1';ones(1,N)];
    p2 = [pt2';ones(1,N)];
    H = p1'\p2';
end

