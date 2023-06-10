%Applies skin_opt.m for mua_der, mus_der, g_der, mua_bone, mus_bone, &
%g_bone
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
lambda=700;             %Wavelength
L = .12;                %Thickness of first layer
a = 0;
b = 11800;
n = 5e5;
ds = (b-a)./n;
s = (a:ds:b);
%------------vvv MAY NEED TO CHANGE BACK
mus_der = 90;
mua_der = 0.02;
mus_bone = 110;
mua_bone = 0.1;
g_der = 0.9;
g_bone = g_der;

Z = 1e-9*(4.62141613 - 0.00034*lambda);

z = 0;

X1 =(142260.46875 - 1523.922119*lambda + 4.152612*lambda.^2)./...
    (1 + 514.571228*lambda - 5.243055*lambda.^2 + 0.013448*lambda.^3 ...
       + 0.000001*lambda.^4 - Z.*lambda.^5);

%mua_der=10.^X1;        %Absorption Coefficient of Dermis

X2 = 1.95108 + .013219*lambda - 4.75777e-5 *lambda.^2 + ...
     6.280311e-8*lambda.^3 - 2.8417e-11*lambda.^4;

%mus_der=10.^X2;        %Scattering Coefficient of Dermis

%g_der=0.81*ones(size(lambda));      %Anisotropy of Dermis

%-------vvv CHANGE TO BONE VALUES 
%mua_bone=427.36./(1+0.059858*(lambda-300))+20;   %Absorption Coeff of Epidermis

X3 =33.1275 - 0.343538*lambda + 1.58018e-3*lambda.^2 -...
    3.78285e-6*lambda.^3 +  4.97727E-9*lambda.^4 -...
    3.42073e-12*lambda.^5 + 9.61503e-16*lambda.^6;

%mus_bone = 10.^X3;                     %Scattering Coefficient of Epidermis
%g_bone = 2.9166e-4*lambda + 0.6225;  %Anisotropy of Epidermis

%--------------------------------------
%EQUIVALENT TRANSPORT OR REDUCED SCATTERING COEFFICIENT
mus_derprime = mus_der.*(1-g_der);
mus_boneprime = mus_bone.*(1-g_bone);
%DIFFUSION CONSTANT
D_der = (3.*(mus_derprime+mua_der))^-1;
D_bone = (3.*(mus_boneprime+mua_bone))^-1;

z_b = 2*D_der;

%TRANSPORT INTERACTION COEFFICIENT
mut_derprime = mua_der+mus_derprime;
mut_boneprime = mua_bone+mus_boneprime;
%TRANSPORT SCATTERING ALBEDO
a_derprime = mus_derprime/(mua_der+mus_derprime);
a_boneprime = mus_boneprime/(mua_bone+mus_boneprime);
%EQUATION FOR EQUIVALENT POINT SOURCE
z_1 = (1-((mut_derprime.*L)+1).*exp(-mut_derprime.*L))./(mut_derprime.*(1-exp(-mut_derprime.*L)));
z_2 = L+(1./mut_boneprime);
%EQUATION FOR STRENGTH OF EQUIVALENT POINT SOURCE
w_der = a_derprime.*(1-exp(-mut_derprime*L));
w_bone = a_boneprime.*exp(-mut_derprime*L);

%(alpha_i)^2 = ((D_1*s^2)+mu_ai)/D_i = s^2+mueff_i
%   so alpha_i = sqrt(((D_1*s^2)+mu_ai)/D_i) = sqrt(s^2+mueff_i)
%where mueff_i = sqrt((mu_ai/D_i))
alpha_der = sqrt(((D_der.*(s.^2))+mua_der)./D_der);
alpha_bone = sqrt(((D_bone.*(s.^2))+mua_bone)./D_bone);
if z <= L
    lilphi_der = (w_der./(2.*alpha_der.*D_der)).*((((alpha_bone.*D_bone)-(alpha_der.*D_der))./(((alpha_bone.*D_bone)+(alpha_der.*D_der)).*exp(alpha_der.*(z+(2.*L)+(2.*z_b)-z_1))))+(exp(-alpha_der.*abs(z-z_1)))-(((alpha_der.*D_der)+(alpha_bone.*D_bone))./(((alpha_der.*D_der)+(alpha_bone.*D_bone)).*(exp(alpha_der.*(z+z_1+(2.*z_b))))+(((alpha_der.*D_der)-(alpha_bone.*D_bone)).*exp(alpha_der.*(z+z_1-(2.*L))))))+(((alpha_der.*D_der)-(alpha_bone.*D_bone))./(((alpha_der.*D_der)+(alpha_bone.*D_bone)).*(exp(alpha_der.*((2.*L)-z_1-z)))+(((alpha_der.*D_der)-(alpha_bone.*D_bone)).*exp(-alpha_der.*(z+z_1+(2.*z_b))))))+(((alpha_bone.*D_bone)-(alpha_der.*D_der))./((((alpha_der.*D_der)+(alpha_bone.*D_bone)).*exp(alpha_der.*(z_1+(2.*z_b)+(2.*L)-z)))+(((alpha_der.*D_der)-(alpha_bone.*D_bone)).*exp(alpha_der.*(z_1-z))))))+((w_bone.*exp(-alpha_bone.*(z_2-L)).*sinh(alpha_der.*(z_b+z)))./((alpha_bone.*D_bone.*sinh(alpha_bone.*(z_b+L)))+(alpha_der.*D_der.*cosh(alpha_der.*(z_b+L)))));

    %lilphi_bone = (((alpha_bone.*D_bone.*w_der.*exp(-alpha_der.*(L-z_1)))-(alpha_der.*D_der.*w_bone.*exp(-alpha_der.*(z_2-L))))./(2.*(alpha_bone.*D_bone.*tanh(alpha_der.*(z_b+L))+(alpha_der.*D_der))))+(((w_bone.*exp(-alpha_bone.*(z_2-L))+(w_der.*exp(-alpha_der.*(L-z_1)))).*sinh(alpha_der.*(L+z_1))-(w_der.*exp(-alpha_der.*(z_b+z_1))))./(2.*((alpha_bone.*D_bone.*sinh(alpha_der.*(z_b+L)))+(alpha_der.*D_der.*cosh(alpha_der.*(z_b+L)))))).*exp(-alpha_bone.*(z-L))+((w_bone.*exp(-alpha_bone.*abs(z_2-z)))./(2.*alpha_bone.*D_bone));
end

dP=0.01;
P=(1:1:100)*dP;
ds=0.1;

value=s.*lilphi_der;            %See Eq. 14
value(value~=value)=0;
%Calculating value of PHI at z = 0 in Eq. 14
for ind= 1:100
     bessel = besselj (0, (s.*P(ind)));
     value = value.*bessel;
     Phi1 = sum(value.*ds);
     PHI1(ind) = Phi1*(1/(2*pi));
end

%Calculates the value of PHI at z = 0 + dz
for k= 1:100
     bessel2 = besselj (0.001, (s.*P(k)));
     value2 = value.*bessel2;
     Phi2 = sum(value2.*ds);
     PHI2(k) = Phi2*(1/(2*pi));
end
 
%Calculates dPHI/dz for Eq 15
%Works as a derivative approximation as (change in PHI)/(change in z)
dPHI = (PHI2-PHI1)./.001;

%EQUATION FOR DIFFUSE REFLECTANCE
R = ((1/4).*PHI1)+((1/2).*D_der.*dPHI);

semilogy (P,R);
 
% rho = 2;
% ds = (0:.01:1e5);
% phi_der = (1/(2.*pi)).*sum(s.*lilphi_der.*besselj(0,srho).*ds);   %Fluence Rate
% 
% for i = 1:2
%     %2 layers: Dermis, Bone
%     
%     %EQUATION FOR DIFFUSIVE REFLECTANCE:
%     R = sum(.25*phi+(.5.*D_der.*dphi./dz));
%     %^^^MUST HAVE VALUES CHANGE FOR APPROPRIATE LAYER
%     display(R)
% end

%QUESTIONS:
%