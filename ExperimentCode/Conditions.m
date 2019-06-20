function [EventTable]=Conditions(reps,cond1,cond2)
i=1;
for x =1:cond1
    for y=1:cond2
        conditions(i,:)=[x y];
        i=i+1;
    end
end
EventTable=[];
for i=1:reps
    EventTable=[EventTable;conditions];
end
Totaltrials=length(EventTable(:,1));
TrialNumber=(1:Totaltrials)';
EventTable=[TrialNumber,EventTable];

randomsequence=randperm(Totaltrials)';
EventTable(:,1)=randomsequence;
EventTable=sortrows(EventTable,1);
end