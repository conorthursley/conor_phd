 %% Function that changes the frequency to find the amplitude ratio 
 %--------------------
 % Purpose of this code is to prepare for f,m, and k optimisation
 % x/X=k2/(k2-m2*w^2)
 tic
 control = 0:1:10;
 hold on
 for ii=1:length(control)
     ratio=control(ii);
     w=linspace(0,25,1000);
        k1=37;
        m1=1;
        k2=ratio*k1;
        m2=0.4647;
        
        disRatio=abs(k2./(k2-m2*w.^2));
        X=disRatio/1000;
        plot(w,X)  
     
 end
 xlabel('Frequency, Hz')
 ylabel('Displacement ratio, m_2/m_1')
 hold off
 grid
 toc
 

 
