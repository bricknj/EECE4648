%
% skin_opt.m
%
% Skin Optical Properties; Unpublished (from memo 1504.)
%
%  Original fits by Jeremy Warner in 1994
%
%
%  Based on data from Anderson and Parish.
%
% Citations:
%
%Van Gemert, M.J.C., Jacques, Steven L., Sterenborg, H.J.C.M.,
%  Star, W.M., ``Skin Optics,'' {\sl IEEE Transactions on Biomedical
%  Engineering}, Dec 1989 v 36 n 12, page 1146.
%  
%  Wan, San, R. Rox Anderson, and John A. Parrish,
%  ``Analytical Modeling for the Optical Properties of the Skin with
%  {\sl In Vitro} and {\sl In Vivo} Applications ,''
%  {\sl Photochemistry and Photobiology 34},
%  Pp.~493--499.
%  
% Especially, see the following one.
%
%  Anderson, R. R., and
%  J. A. Parrish,
%  ``Optical Properties of Human Skin,''
%  in J. D. Regen, and J. A. Parrish, eds.,
%  {\sl The Science of Photomedicine},New York, Plenum, 1982.
%  Pp. 147--194.

%
% Chuck DiMarzio, Northeastern University, Aug 2008
%
%  See also skin_sj.m for a different result
%
lambda=400:10:800;
mua_epi=427.36./(1+0.059858*(lambda-300))+20;

X1 =33.1275 - 0.343538*lambda + 1.58018e-3*lambda.^2 -...
    3.78285e-6*lambda.^3 +  4.97727E-9*lambda.^4 -...
    3.42073e-12*lambda.^5 + 9.61503e-16*lambda.^6;

mus_epi=10.^X1;
g_epi = 2.9166e-4*lambda + 0.6225;

Z = 1e-9*(4.62141613 - 0.00034*lambda);

X2 =(142260.46875 - 1523.922119*lambda + 4.152612*lambda.^2)./...
    (1 + 514.571228*lambda - 5.243055*lambda.^2 + 0.013448*lambda.^3 ...
       + 0.000001*lambda.^4 - Z.*lambda.^5);

mua_der=10.^X2;

X3 = 1.95108 + .013219*lambda - 4.75777e-5 *lambda.^2 + ...
     6.280311e-8*lambda.^3 - 2.8417e-11*lambda.^4;

mus_der=10.^X3;

g_der=0.81*ones(size(lambda));

figure;
subplot(2,1,1);
semilogy(lambda(1:10:end),mua_epi(1:10:end),...
       'b-+',lambda(1:10:end),mus_epi(1:10:end),'c-o',...
       lambda(1:10:end),mua_der(1:10:end),'r-^',...
       lambda(1:10:end),mus_der(1:10:end),'m-v',...
       lambda,mua_epi,'b-',...
       lambda,mus_epi,'c-',...
       lambda,mua_der,'r-',...
       lambda,mus_der,'m-');
legend('\mu_a(e)','\mu_s(e)','\mu_a(d)','\mu_s(d)')
subplot(2,1,2);
plot(lambda(1:10:end),g_epi(1:10:end),'b-+',...
     lambda(1:10:end),g_der(1:10:end),'r-v',...
     lambda,g_epi,'b-',...
     lambda,g_der,'r-');
legend('g(e)','g(d)');


