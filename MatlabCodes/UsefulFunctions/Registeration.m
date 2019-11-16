Session=6:12;
loadPath='Q:\L365\Mapping\20170330\AreaMask';
AddressSetup;
t=Tiff(strcat(tiffIn{Session(1)},'00001.tif'),'r');
FirstFrame=t.read();
 %FirstFrame=imrotate(FirstFrame,90);
t.close();

[optimizer, metric] = imregconfig('multimodal');

optimizer.InitialRadius = 0.001;
optimizer.Epsilon = 1.5e-6;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 1000;
maskRegistered=zeros(1017,1017,8); 
for s=9:10%Session
    
t=Tiff(strcat(tiffIn{s},'00001.tif'),'r');
fixed=t.read();
if (s==13) fixed=imrotate(fixed,90); end
t.close();
    
    
tform = imregtform(FirstFrame, fixed, 'similarity', optimizer, metric);

for i=1:8
    [~,~,mask]=imread(strcat(loadPath,'\a',int2str(i),'.png'));
    mask=(mask ~=0);
    maskRegistered(:,:,i) = imwarp(mask,tform,'OutputView',imref2d(size(mask)));
    imwrite(maskRegistered(:,:,i),strcat(LoadPath{s+20},'\areaMasks\','a',int2str(i),'.png'));
end
end

%%%%% reload after correction and save as matlab array
Session=27:30
for s=Session
    for i=1:8
        x=imread(strcat(LoadPath{s+20},'\areaMasks\a',int2str(i),'.png'));
        Area(:,:,i)=x;%(4:1020,4:1020);
    end
    save(strcat(LoadPath{s+20},'\areas'),'Area');
end







