%Generate uniform integers on the set 1:n:

N=204;%94
n=195;
filename=sprintf('group-%d.txt',n);
%p=randperm(N);
r=zeros(1,n);
r(1)=rand;
r(1)=ceil(N*r(1));
i=2;
while (i<=n)
    r(i)=rand;
    r(i)=ceil(N*r(i));
    for j=1:i-1
        if (r(i)==r(j)) 
            i=i-1;
            break;
        end
    end
    i=i+1;
end
r=sort(r);
dlmwrite(filename,r,' ');
