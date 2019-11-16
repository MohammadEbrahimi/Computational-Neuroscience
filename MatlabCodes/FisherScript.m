
LoadPath{17}= 'C:\Users\Sadegh\Documents\VLMReborn\L347\Data\allDays';
LoadPath{18}= 'C:\Users\Sadegh\Documents\VLMReborn\L354\Data\allDays';
LoadPath{19}= 'C:\Users\Sadegh\Documents\VLMReborn\L362bis\Data\allDays';
LoadPath{20}= 'C:\Users\Sadegh\Documents\VLMReborn\L364\Data\allDays';
LoadPath{21}= 'C:\Users\Sadegh\Documents\VLMReborn\L365\Data\allDays';
LoadPath{22}= 'C:\Users\Sadegh\Documents\VLMReborn\L367\Data\allDays';
LoadPath{23}= 'C:\Users\Sadegh\Documents\VLMReborn\L368\Data\allDays';



train=1;

Days=[17:23];
%%%%%
result=zeros(2,7);
for i=1:length(Days)
    Day=Days(i)
    result(1,i)=DatasetSize(LoadPath{Day},'M','H',1);
    result(2,i)=DatasetSize(LoadPath{Day},'M','H',2);
end







