function [den,freq]= density_function(x)

% den=zeros(1,51);
% freq=zeros(1,4);
% num=zeros(1,4); % A,C,G,U
% 
% for i=1:51
%     k=4*(i-1);
%     if x(k+1)==1&&x(k+2)==0&&x(k+3)==0&&x(k+4)==0
%         num(1)= num(1)+1;
%         den(i)=num(1)/i;
%     end
%     
%     if x(k+1)==0&&x(k+2)==1&&x(k+3)==0&&x(k+4)==0
%         num(2)= num(2)+1;
%         den(i)=num(2)/i;
%     end
%     
%     if x(k+1)==0&&x(k+2)==0&&x(k+3)==1&&x(k+4)==0
%         num(3)= num(3)+1;
%         den(i)=num(3)/i;
%     end
%     
%     if x(k+1)==0&&x(k+2)==0&&x(k+3)==0&&x(k+4)==1
%         num(4)= num(4)+1;
%         den(i)=num(4)/i;
%     end
% end
% 
% for i=1:4
%     freq(i)= num(i)/51;
% end
%     
% return

den=zeros(1,51);
freq=zeros(1,4);
num=zeros(1,4); % A,C,G,U

for i=1:51
    k=3*(i-1);
    if x(k+1)==1&&x(k+2)==1&&x(k+3)==1
        num(1)= num(1)+1;
        den(i)=num(1)/i;
    end
    
    if x(k+1)==0&&x(k+2)==0&&x(k+3)==1
        num(2)= num(2)+1;
        den(i)=num(2)/i;
    end
    
    if x(k+1)==1&&x(k+2)==0&&x(k+3)==0
        num(3)= num(3)+1;
        den(i)=num(3)/i;
    end
    
    if x(k+1)==0&&x(k+2)==1&&x(k+3)==0
        num(4)= num(4)+1;
        den(i)=num(4)/i;
    end
end

for i=1:4
    freq(i)= num(i)/51;
end
    
return

