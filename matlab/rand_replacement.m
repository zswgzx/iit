%create random num from a line in txt file with replacement
clear all;close all;clc;

n=195;N=100;      %n subjects, N combinations
flag=0;         %determine if it is repeated
file=sprintf('group-%d.txt',n);
a=importdata(file);
list=zeros(N,n);

for i=1:N
    random=ceil(n.*rand(1,n));
    random=sort(random);
    list(i,:)=random;
    if i>1
     for j=1:i-1
           compare=list(i,:)-list(j,:);
          if (nnz(compare)==0) 
              %duplication detected
            i=i-1;
            flag=1;
          else
            flag=0;
          end
     end
    end
    
    if flag==0
        for j=1:n
           list(i,j)=a(random(j));
        end
    end
end

%write to a txt file
file=sprintf('%d_100rep.txt',n);
fp_txt=fopen(file,'w');
    for i=1:N
        for j=1:n
            fprintf(fp_txt,'%d',list(i,j));
            if (j==n)
                fprintf(fp_txt,'\n');
            else
                fprintf(fp_txt,' ');
            end
        end
    end
fclose(fp_txt);
