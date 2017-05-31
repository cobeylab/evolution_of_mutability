function p_m_mer=prob_Lmer(h, f_0_mer)
%calculates the probabilty of each L-mer according to background frequency
%f_0_mer (vector length 4^L) and model scores h (matrix L*4).
L=7;
%%
ii=zeros(4^L,L); %matrix that holds the indentity of the bases in each position (L) for each L-mer (L^4).
for i=1:length(f_0_mer(:))
    t=str2num(dec2base(i-1,4,L)')+1; %convert the number of the L-mer to the actual L-mer content
    ii(i,:)=t(end:-1:1); %reverse
end
%%
E=zeros(size(f_0_mer)); %score for each L-mer to be calculated
for i=1:length(f_0_mer(:))
    E(i)=0;
    for f=1:L
        E(i)=E(i)+h(f,ii(i,f)); %adding score acording to model h and indentity of bases for current L-mer
    end
end
p_m_mer=f_0_mer.*exp(E); %multiply score according to model and background L-mer frequency
p_m_mer=p_m_mer/sum(p_m_mer(:)); %normalize probabilities across L-mer