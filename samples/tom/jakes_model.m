close all; clear all;

% fm Doppler frequencies is Hz
fm=[1 10 100]; 

M=7; %seven scatterers
t=0:0.001:2; %0 to 2 sec, 1ms sample period

h = zeros( length(fm), length(t) );

for i=1:3
    h = jakes_response(fm(i), t, M);
    
    r=abs(h);
    
    %find the envelope average
    mean=sum(r)/length(t);
    subplot(3,1,i);
    
    %plot the figure and shift the envelope average to 0dB
    plot(t,(10*log10(r)-10*log10(mean)));
    titlename=['fd = ' int2str(fm(i)) ' Hz'];
    title(titlename);
    xlabel('time (s)');
    ylabel('Gain (dB)');
end

function h = jakes_response(fm, t, M)
% fm - doppler frequency of channel in Hz
% t - vector indicating time location of samples in sec)
% M - number of scatterers (7 is typically sufficient)

    Wm=2*pi*fm; % angular freq
    N = 4*M+2; % defines angular spacing of scatterers
    Wn(M)=0; beta(M)=0;
    
    ritemp(M,length(t))=0; 
    rqtemp(M,length(t))=0; 
    rialpha(1,length(t))=0; 
    
    for n=1:1:M
        gamma = pi*rand(1,M); %each scatterer has a random phase
        
        for i=1:length(t)
            Wn(n)=Wm*cos(2*pi*n/N);
            beta(n)=pi*n/M;
            ritemp(n,i)=2*cos(beta(n))*cos(Wn(n)*t(i) + gamma(n));
            rqtemp(n,i)=2*sin(beta(n))*cos(Wn(n)*t(i) + gamma(n));
            rialpha(1,i)=2*cos(Wm*t(i))/sqrt(2);
        end
    end
    
    %Add scatterers together coherently
    ri=sum(ritemp)+rialpha;
    rq=sum(rqtemp);
    h = ri + 1i*rq;
end