clear all;
img = load_img();
imgout = img_stitch(img{1}, img{2});
imgout = img_stitch(imgout, img{3});
imgout = img_stitch(imgout, img{4});
imshow(imgout);
imwrite(imgout,'out.jpg','jpg');
