function img = load_img()
    disp('Loading images ...');
    img_list = dir('./input_s/C*.jpg');
    img = cell(1,length(img_list));
    for i = 1:length(img_list) 
        img{i} = imread(['./input_s/',img_list(i).name]);
    end
    disp('images loaded!');
end

