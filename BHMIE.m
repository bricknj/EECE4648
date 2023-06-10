function [QEXT,QSCA,QBACK,QSCATG,G,S1,S2]=BHMIE(X,REFREL,NANG)
%
% bhmie.m
%  by Tom Gaudette at Northeastern University, late 1990's
%
%  Boren and Huffman Mie Scattering code.
%  [QEXT,QSCA,QBACK,QSCATG,G,S1,S2]=BHMIE(X,REFREL,NANG)
% 
%  x is 2*pi*radius*index_of_background/wavelength = kr, wavelength is in the 
%       background medium
%  refrel is the ratio of the particle index to that of the background
%  nang is the number of angles for which to plot the results
%  Qext is the extinction efficiency,
%  Qsca is the scattering efficiency,
%  Qback is the backscattering efficiency,
%  Qscatg is ...
%  G is the anisotropy
%  S1 and S2 are scattering amplitudes; see Boren and Huffman
%
%  To get cross sections, multiply Q by \pi radius^2
%
%  Modified by Jose Barbosa 21 July 2002, corrected comment on x 
%           parameter
%
% $Version 1.0 
DX=X;
Y= X*REFREL;
XSTOP=X+4.*X^0.33333+2.0;
NSTOP=XSTOP;
YMOD=abs(Y);
NMX=ceil(max(XSTOP,YMOD)+15);
DANG=(pi/2)/(NANG-1);
      
%      for j=1:NANG;
%      	THETA(J)=(J-1.0)*DANG;
%    	   AMU(J)=COS(THETA(J));
%  	 end	
theta=[0:1:(NANG-1)]*DANG;
AMU=cos(theta);
%
%Logarithmic derivative D(J) calculated by downward
%recurrence beginning with initial value 0.0 +1*0.0 at j=nmx
%
D=zeros(1,NMX);
NN=NMX-1;
for N=1:NN;
   RN=NMX-N+1;
   % 120  D(NMX-N)=(RN/Y)-(1./(D(RN)+RN/Y));
   D(NMX-N)=(RN/Y)-(1./(D(NMX-N+1)+RN/Y));
end
%      for J=1:NANG;
%      	PI1(J)=1.0;
%		   PI0(J)=0.0;
%	    end
PI1=ones(NANG,1);
PI0=zeros(NANG,1);
	
%      NN=2*NANG-1
%      for J=1:NN
% 		   S1(J)=CMPLX(0.0,0.0);
% 		   S2(J)=CMPLX(0.0,0.0);
%	    end
S1=zeros(2*NANG-1,1)+i*zeros(2*NANG-1,1);
S2=zeros(2*NANG-1,1)+i*zeros(2*NANG-1,1);

%
%     Riccati-bessel functions with real argument x calculated
%     by upward recurrence
%
%i=sqrt(-1);
PSI0=cos(DX);
PSI1=sin(DX);
CHI0=-sin(X);
CHI1=cos(X);
APSI0=PSI0;
APSI1=PSI1;
XI0=APSI0-i*CHI0;
XI1=APSI1-i*CHI1;
QSCA=0.0;
QSCATG=0.0;
N=1;
ANOLD=0;
BNOLD=0;

Tom=1;     
while(Tom)
   DN=N;
   RN=N;
   FN= (2*RN+1)/(RN*(RN+1.0));
   PSI=(2*DN-1)*PSI1/DX-PSI0;
   APSI=PSI;
   CHI=(2*RN-1)*CHI1/X-CHI0;
   XI =APSI-i*CHI;
   AN=(D(N)/REFREL+RN/X)*APSI-APSI1;
   AN=AN/((D(N)/REFREL+RN/X)*XI-XI1);
   BN=(D(N)*REFREL+RN/X)*APSI-APSI1;
   BN=BN/((D(N)*REFREL+RN/X)*XI-XI1);
   QSCA=QSCA+(2.*RN+1.0)*(abs(AN)^2+abs(BN)^2);
   for J=1:NANG;
     	JJ=2*NANG-J;
     	PI(J)=PI1(J);
      TAU(J)=RN*AMU(J)*PI(J)-(RN+1.0)*PI0(J);
      P=(-1.0)^(N-1);
      S1(J)=S1(J)+FN*(AN*PI(J) +BN*TAU(J));
      T=(-1.0)^N;
      S2(J)=S2(J)+FN*(AN*TAU(J) +BN*PI(J));
      if(J~=JJ)
           S1(JJ)=S1(JJ)+FN*(AN*PI(J)*P+BN*TAU(J)*T);
           S2(JJ)=S2(JJ)+FN*(AN*TAU(J)*T +BN*PI(J)*P);
      end;
    end;
	  
    AOASTR = ANOLD * conj(AN);
    BOBSTR = BNOLD * conj(BN);
    ABSTAR = AN    * conj(BN);
    QSCATG = QSCATG + ((RN-1.)*(RN+1.)/(RN)) * real(AOASTR+BOBSTR) + FN * real(ABSTAR);
    
    ANOLD = AN;
    BNOLD = BN;
    PSI0=PSI1;
    PSI1=PSI;
    CHI0=CHI1;
    CHI1=CHI;
    APSI1=PSI1;
    XI1=APSI1-i*CHI1;
    N=1+N;
    RN=N;
	  
    for J=1:NANG;
      	PI1(J)=((2.*RN- 1.)/(RN-1.0))*AMU(J)*PI(J);
         PI1(J)=PI1(J)-RN*PI0(J)/(RN-1.0);
         PI0(J)=PI(J);
     end;
   %    IF(N-1-NSTOP) 200,300,300
    A_tom=N-1-NSTOP;
    if (A_tom>=0)
       Tom=0;
       break;
    end;
end;
      
      QEXT=(4.0/X^2)* real(S1(1));
      QSCA = (2.0/X^2)*QSCA;
      QBACK=(4.0/X^2)*abs(S1(2*NANG-1))*abs(S1(2*NANG-1));
      QSCATG = (4.0/X^2)*QSCATG;
      G = QSCATG/QSCA;
