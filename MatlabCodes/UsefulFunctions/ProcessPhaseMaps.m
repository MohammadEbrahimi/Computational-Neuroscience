Mouse='L364';
AddressSetup;
if Mouse=='L364'
    add1='Q:\L364\Mapping\20170328mapping\';
    add2='Q:\L364\Mapping\20170328mapping\';
    firstSession=1;
    addff='G:\2017-03-11\L364\Calcium\trial1\L3641_00001.tif';
elseif Mouse=='L365'
    add1='Q:\L365\Mapping\20170330\';
    add2='Q:\L365\Mapping\20170411\';
    firstSession=6;
    addff='G:\2017-03-09\L365\Calcium\trial1\L3651_00001.tif';
elseif Mouse=='L367'
    add1='Q:\L367\Mapping\20170419-21\';
    add2='Q:\L367\Mapping\20170419-21\';
    firstSession=13;
    addff='G:\2017-03-06\L367\Calcium\trial1\L3671_00001.tif';
elseif Mouse=='L368'
    add1='Q:\L368\Mapping\20170330\';
    add2='Q:\L368\Mapping\20170330\';
    firstSession=20;
    addff='G:\2017-03-09\L368\Calcium\trial1\L3681_00001.tif';
end


D0=open(strcat(add1,'Afp0.mat'));
D90=open(strcat(add1,'Afp90.mat'));
D180=open(strcat(add2,'Afp180.mat'));
D270=open(strcat(add2,'Afp270.mat'));
t=Tiff(strcat(add1,Mouse,'deg01_00001.tif'),'r');
ff0=t.read();
t.close();
t=Tiff(strcat(add1,Mouse,'deg901_00001.tif'),'r');
ff90=t.read();
t.close();
t=Tiff(strcat(add2,Mouse,'deg1801_00001.tif'),'r');
ff180=t.read();
t.close();
t=Tiff(strcat(add2,Mouse,'deg2701_00001.tif'),'r');
ff270=t.read();
t.close();
t=Tiff(addff,'r');
FirstFrame=t.read();
% FirstFrame=imrotate(FirstFrame,90);
t.close();



    
    [optimizer, metric] = imregconfig('multimodal');
    Dim=[1017,1017];
    optimizer.InitialRadius = 0.001;
    optimizer.Epsilon = 1.5e-6;
    optimizer.GrowthFactor = 1.01;
    optimizer.MaximumIterations = 1000;
    
    tform0 = imregtform(ff0,FirstFrame, 'similarity', optimizer, metric);
   tform90 = imregtform(ff90, FirstFrame, 'similarity', optimizer, metric);
    tform180 = imregtform(ff180,FirstFrame, 'similarity', optimizer, metric);
    tform270 = imregtform(ff270,FirstFrame, 'similarity', optimizer, metric);

    p0 = imwarp(imresize(D0.p,[1024,1024]),tform0,'OutputView',imref2d(Dim));
       p90 = imwarp(imresize(D90.p,[1024,1024]),tform0,'OutputView',imref2d(Dim));
    p180 = imwarp(imresize(D180.p,[1024,1024]),tform0,'OutputView',imref2d(Dim));
    p270 = imwarp(imresize(D270.p,[1024,1024]),tform0,'OutputView',imref2d(Dim));

    A0 = imwarp(imresize(D0.A,[1024,1024]),tform0,'OutputView',imref2d(Dim));
     A90 = imwarp(imresize(D90.A,[1024,1024]),tform0,'OutputView',imref2d(Dim));
    A180 = imwarp(imresize(D180.A,[1024,1024]),tform0,'OutputView',imref2d(Dim));
    A270 = imwarp(imresize(D270.A,[1024,1024]),tform0,'OutputView',imref2d(Dim));


         figure
         ffreg=imwarp(ff0,tform0,'OutputView',imref2d(Dim));
         imshowpair(FirstFrame, ffreg,'Scaling','joint');
         

A=A0+A90+A180+A270;

imwrite(ind2rgb(p0,jet(44)),strcat(add1,'p0.png'));
imwrite(ind2rgb(p90,jet(34)),strcat(add1,'p90.png'));
imwrite(ind2rgb(p180,jet(44)),strcat(add2,'p180.png'));
imwrite(ind2rgb(p270,jet(83)),strcat(add2,'p270.png'));
imwrite(A*10,strcat(add1,'A.png'));




