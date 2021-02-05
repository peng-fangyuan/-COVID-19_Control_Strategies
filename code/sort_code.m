%get the parameter
para=squeeze(median(paraall(:,:,:,10),2));
pa=zeros(40,1);
pa(1:8)=mean(para(:,1:14),2);
pa(9:16)=mean(para(:,15:24),2);
pa(17:24)=mean(para(:,25:34),2);
pa(25:32)=mean(para(:,35:44),2);
pa(33:40)=mean(para(:,45:end),2);