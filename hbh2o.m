% transmisson of blood and water
% Chuck DiMarzio, Northeastern University, July 2004
w=[250:2:1000];
[ko,kd]=hbspec(w);
ko=ko/4/log(10);;
kd=kd/4/log(10);;
for qqq=1:length(w);
[ri(qqq),ii(qqq)]=indwat(w(qqq)/1000,0.);
end;
ext=2*pi*ii./(w*1e-7);
figure;semilogy(w,ko*10*1.2e-4,'r-',w,kd*10*1.2e-4,'b-',...
    w,ko*10*23e-4,'m-',w,kd*10*23e-4,'c-',...
    w,ext,'g-',...
    w(1:10:end),ko(1:10:end)*10*1.2e-4+ext(1:10:end),'ko')
grid on;
xlabel('\lambda,Wavelength, nm');ylabel('\mu_a, Absorption Coeff, /cm');
legend('Oxy in Skin','Deoxy','Oxy in Blood','Deoxy','Water');



