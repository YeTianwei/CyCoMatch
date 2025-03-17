function [k] = findK(D)

    d = D/mean(D(1:100));
    k=length(d);
    for i=1:length(d)-10
        if sum(abs(d(i:i+10)))<0.1
            k = i;
            break;
        end
    end
end