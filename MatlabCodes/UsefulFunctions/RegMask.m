tiffIn{1}='G:\2017-03-11\L364\Calcium\trial1\L3641_';
tiffIn{2}='G:\2017-03-12\L364\Calcium\L364\L3641_';
tiffIn{3}='G:\2017-03-14\L364\Calcium\20170314\L3641_';
tiffIn{4}='G:\2017-03-16\L364\Calcium\trial3&4\L364part4\L3641_';
tiffIn{5}='G:\2017-03-17\L364\Calcium\20170317\L3641_';
tiffIn{6}='G:\2017-03-09\L365\Calcium\trial1\L3651_';
tiffIn{7}='G:\2017-03-10\L365\Calcium\trial1\L3651_';
tiffIn{8}='G:\2017-03-11\L365\Calcium\trial1\L3651_';
tiffIn{9}='G:\2017-03-12\L365\Calcium\L365\L3651_';
tiffIn{10}='G:\2017-03-12\L365\Calcium\L365\L365part21_';
tiffIn{11}='G:\2017-03-16\L365\Calcium\20170316\L3651_';
tiffIn{12}='G:\2017-03-17\L365\Calcium\part2\L36521_';
tiffIn{13}='G:\2017-03-06\L367\Calcium\trial1\L3671_';
tiffIn{14}='G:\2017-03-07\L367\Calcium\trial1\L3671_';
tiffIn{15}='G:\2017-03-08\L367\Calcium\trial1\L3671_';
tiffIn{16}='G:\2017-03-09\L367\Calcium\trial1\L3671_';
tiffIn{17}='G:\2017-03-09\L367\Calcium\trial1\L367p21_';
tiffIn{18}='G:\2017-03-11\L367\Calcium\trial1\L3671_';
tiffIn{19}='G:\2017-03-11\L367\Calcium\trial1\L367part21_';
tiffIn{20}='G:\2017-03-09\L368\Calcium\trial1\L3681_';
tiffIn{21}='G:\2017-03-10\L368\Calcium\trial1\L3681_';
tiffIn{22}='G:\2017-03-11\L368\Calcium\trial1\L3682_';
tiffIn{23}='G:\2017-03-12\L368\Calcium\20170312\L3681_';
tiffIn{24}='G:\2017-03-14\L368\Calcium\20170314\L3681_';
tiffIn{25}='G:\2017-03-16\L368\Calcium\20170316\L3681_';
tiffIn{26}='G:\2017-03-17\L368\Calcium\20170317\2\L368p21_';


CaMovies{1}='Q:\L364\Calcium\2017_03_11\';
CaMovies{2}='Q:\L364\Calcium\2017_03_12\';
CaMovies{3}='Q:\L364\Calcium\2017_03_14\';
CaMovies{4}='Q:\L364\Calcium\2017_03_16\';
CaMovies{5}='Q:\L364\Calcium\2017_03_17\';
CaMovies{6}='Q:\L365\Calcium\2017_03_09\';
CaMovies{7}='Q:\L365\Calcium\2017_03_10\';
CaMovies{8}='Q:\L365\Calcium\2017_03_11\';
CaMovies{9}='Q:\L365\Calcium\2017_03_12_s1\';
CaMovies{10}='Q:\L365\Calcium\2017_03_12_s2\';
CaMovies{11}='Q:\L365\Calcium\2017_03_16\';
CaMovies{12}='Q:\L365\Calcium\2017_03_17\';
CaMovies{13}='Q:\L367\Calcium\2017_03_06\';
CaMovies{14}='Q:\L367\Calcium\2017_03_07\';
CaMovies{15}='Q:\L367\Calcium\2017_03_08\';
CaMovies{16}='Q:\L367\Calcium\2017_03_09_s1\';
CaMovies{17}='Q:\L367\Calcium\2017_03_09_s2\';
CaMovies{18}='Q:\L367\Calcium\2017_03_11_s1\';
CaMovies{19}='Q:\L367\Calcium\2017_03_11_s2\';
CaMovies{20}='Q:\L368\Calcium\2017_03_09\';
CaMovies{21}='Q:\L368\Calcium\2017_03_10\';
CaMovies{22}='Q:\L368\Calcium\2017_03_11\';
CaMovies{23}='Q:\L368\Calcium\2017_03_12\';
CaMovies{24}='Q:\L368\Calcium\2017_03_14\';
CaMovies{25}='Q:\L368\Calcium\2017_03_16\';
CaMovies{26}='Q:\L368\Calcium\2017_03_17\';


Sessions=20:26;
t=Tiff(strcat(tiffIn{Sessions(1)},num2str(1,'%05d'),'.tif'),'r');
FirstFrameRef=t.read();
t.close();


[~,~,Mask]=imread(strcat(CaMovies{Sessions(1)},'L368_Mask.png'));

for s=Sessions
    

t=Tiff(strcat(tiffIn{s},num2str(1,'%05d'),'.tif'),'r');
FirstFrame=t.read();
t.close();



[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = 0.001;
optimizer.Epsilon = 1.5e-6;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 1000;

tform = imregtform(FirstFrameRef, FirstFrame, 'similarity', optimizer, metric);

RegisteredMask = imwarp(Mask,tform,'OutputView',imref2d(size(FirstFrame)));
imwrite((RegisteredMask~=0),strcat(CaMovies{s},'S',int2str(s),'_Mask.png'));
RegisteredMask= uint16(RegisteredMask~=0) .*(max(max(FirstFrame))*2);
figure();

imshowpair(RegisteredMask, FirstFrame,'Scaling','joint');
title(strcat('Session',num2str(s)));
end


