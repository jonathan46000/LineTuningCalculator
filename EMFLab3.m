%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    EMF Lab 3: Impedance Matching
%    Jodlet Saintcyr and Jonathan Lundquist
%    March 29, 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%_INPUTS_%%%

prompt = 'Please Enter Characteristic Impedance of Tx Line (Number):';
Zo = input(prompt);

prompt = 'Please Enter the Real Component of Load Impedance (Number):';
ZL = input(prompt);

prompt = 'Please Enter the Imaginary Component of Load Impedance (Number):';
ZL = ZL + i*input(prompt);

prompt = 'Please Enter the Frequency of Operation (A value between 1 and 6 GHz):';
f = input(prompt)*(1e9);
    while (f < (1e9)) || (f > (6e9))
        prompt = 'Please enter a value in GHz between 1 and 6:';
        f = input(prompt)*(1e9);
    end
    
prompt = 'Please Enter Relative Permittivity of the Line (Number)';
Er = input(prompt);

%%%_CONSTANTS_%%%

Eo = 8.85e-12;
Uo = 4*pi*(1e-7);

%%%_DERIVED_VALUES_%%%

phase_velocity = 1/sqrt(Er*Eo*Uo);
wavelength = phase_velocity/f;
Yo = 1/Zo;
B = 2*pi/wavelength;
Refl = (ZL - Zo)/(ZL + Zo);
RefMag = abs(Refl);
RefAng = angle(Refl);

%%%_VARIABLES_AND_ARRAYS_%%%

ARRAY_LENGTH = 100;

low_band_freq_1 = 1e9;
low_band_freq_2 = 6e9;

high_band_freq_1 = 24e9;
high_band_freq_2 = 85e9;

low_freq = low_band_freq_1:(floor((low_band_freq_2-low_band_freq_1)/ARRAY_LENGTH)):low_band_freq_2;
high_freq = high_band_freq_1:(floor((high_band_freq_2-high_band_freq_1)/ARRAY_LENGTH)):high_band_freq_2;

%%%_STUB_TUNING_%%%

   %%%_SHUNT_STUB_TUNING_%%%
   
    
      t1 = (imag(ZL) + sqrt(((real(ZL)/(Zo)))*(((real(ZL)-Zo)^2)+(imag(ZL))^2)))/(real(ZL)-Zo);
      t2 = (imag(ZL) - sqrt(((real(ZL)/(Zo)))*(((real(ZL)-Zo)^2)+(imag(ZL))^2)))/(real(ZL)-Zo);
      if isnan(t1)
          t1 = 0;
      end
      if isnan(t2)
          t2 = 0;
      end
      d1 = atan(t1)/B;
      d2 = atan(t2)/B;
      if (d1 < 0)
          d1 = d1 + wavelength/2;
      end
      if (d2 < 0)
          d2 = d2 + wavelength/2;
      end
      
      Bin1 = Yo*(((real(ZL))^2)*t1-(Zo-imag(ZL)*t1)*(Zo*t1+imag(ZL)))/(((real(ZL))^2)+(imag(ZL)+Zo*t1)^2);
      Bin2 = Yo*(((real(ZL))^2)*t2-(Zo-imag(ZL)*t2)*(Zo*t2+imag(ZL)))/(((real(ZL))^2)+(imag(ZL)+Zo*t2)^2);
      
      Bs1 = -Bin1;
      Bs2 = -Bin2;
      
      if(isnan(Bs1))
          Bs1 = 0;
      end
      if(isnan(Bs2))
          Bs2 = 0;
      end
   
      %%%_OPEN_%%%
         
    
      l_stub_open1 = wavelength*(atan(Bs1/Yo))/(2*pi);
      l_stub_open2 = wavelength*(atan(Bs2/Yo))/(2*pi);
             
      if (l_stub_open1 < 0)
          l_stub_open1 = l_stub_open1 + wavelength/2;
      end
      if (l_stub_open2 < 0)
          l_stub_open2 = l_stub_open2 + wavelength/2;
      end

      l_stub_short1 = -wavelength*(atan(Yo/Bs1))/(2*pi);
      l_stub_short2 = -wavelength*(atan(Yo/Bs2))/(2*pi);
             
      if (l_stub_short1 < 0)
         l_stub_short1 = l_stub_short1 + wavelength/2;
      end
      if (l_stub_short2 < 0)
         l_stub_short2 = l_stub_short2 + wavelength/2;
      end

      
   %%%_SERIES_STUB_TUNING_%%%
   
      YL = 1/ZL;
      Yo = 1/Zo;
      
      if (real(YL) == Yo)
         
         t3 = -imag(YL)/(2*Yo);
         t4 = t3;
      else
         t3 = (imag(YL) + sqrt(((real(YL)/(Yo)))*(((real(YL)-Yo)^2)+(imag(YL))^2)))/(real(YL)-Yo);
         t4 = (imag(YL) - sqrt(((real(YL)/(Yo)))*(((real(YL)-Yo)^2)+(imag(YL))^2)))/(real(YL)-Yo);
      end    
      d3 = atan(t3)/B;
      d4 = atan(t4)/B;
      if (d3 < 0)
          d3 = d3 + wavelength/2;
      end
      if (d4 < 0)
          d4 = d4 + wavelength/2;
      end
      
      Xin3 = (((real(YL))^2)*t3-(Yo-t3*imag(YL))*(imag(YL)+t3*Yo))/(Yo*(((real(YL))^2)+(imag(YL)+Yo*t3)^2));
      Xin4 = (((real(YL))^2)*t4-(Yo-t4*imag(YL))*(imag(YL)+t4*Yo))/(Yo*(((real(YL))^2)+(imag(YL)+Yo*t4)^2));
      
      Xs3 = -Xin3;
      Xs4 = -Xin4;
      
      if(isnan(Xs3))
          Xs3 = 0;
      end
      if(isnan(Xs4))
          Xs4 = 0;
      end
   
      %%%_OPEN_%%%
      
      l_series_open1 = -wavelength*(atan(Zo/Xs3))/(2*pi);
      l_series_open2 = -wavelength*(atan(Zo/Xs4))/(2*pi);
             
      if (l_series_open1 < 0)
          l_series_open1 = l_series_open1 + wavelength/2;
      end
      if (l_series_open2 < 0)
          l_series_open2 = l_series_open2 + wavelength/2;
      end

      l_series_short1 = wavelength*(atan(Xs3/Zo))/(2*pi);
      l_series_short2 = wavelength*(atan(Xs4/Zo))/(2*pi);
             
      if (l_series_short1 < 0)
          l_series_short1 = l_series_short1 + wavelength/2;
      end
      if (l_series_short2 < 0)
          l_series_short2 = l_series_short2 + wavelength/2;
      end
      
%%%_QUARTER_WAVE_TRANSFORMER_%%%
   
Zmax = -(2*pi - RefAng)/(2*B);

while (Zmax < 0.000)
    Zmax = Zmax + wavelength/2;
end

if (Zmax > wavelength/4)
    Zmin = Zmax - wavelength/4;
else
    Zmin = Zmax + wavelength/4;
end

ZmaxIn = Zo*((real(ZL) + i*imag(ZL)) + i*Zo*tan(B*Zmax))/(Zo + i*(real(ZL) + i*imag(ZL))*tan(B*Zmax));
ZminIn = Zo*((real(ZL) + i*imag(ZL)) + i*Zo*tan(B*Zmin))/(Zo + i*(real(ZL) + i*imag(ZL))*tan(B*Zmin));

if (abs(imag(ZmaxIn)) < 1e-5)
    ZmaxIn = real(ZmaxIn);
else
    disp("Error Your Math is Borken");
end

if (abs(imag(ZminIn)) < 1e-5)
    ZminIn = real(ZminIn);
else
    disp("Error Your Math is Borken");
end

ZoPrimeMin = sqrt(Zo*ZminIn);
ZoPrimeMax = sqrt(Zo*ZmaxIn);

QuarterWaveZinMaxSolution = ZoPrimeMax*(ZmaxIn + i*ZoPrimeMax*tan(B*wavelength/4))/(ZoPrimeMax+i*ZmaxIn*tan(B*wavelength/4));
QuarterWaveZinMinSolution = ZoPrimeMin*(ZminIn + i*ZoPrimeMin*tan(B*wavelength/4))/(ZoPrimeMin+i*ZminIn*tan(B*wavelength/4));

%%%_DOUBLE_STUB_TUNING_%%%

   % To be implemented at some future date

%%%_IMPEDANCE&REFLECTION_AS_A_FUNCTION_OF_FREQ_%%%

   %%%_Z(f)_DOWNLINE_OF_TUNERS_%%%
   
   Zd1 = input_impedance(Zo,ZL,phase_velocity./low_freq,d1);
   Zd2 = input_impedance(Zo,ZL,phase_velocity./low_freq,d2);
   Zd3 = input_impedance(Zo,ZL,phase_velocity./low_freq,d3);
   Zd4 = input_impedance(Zo,ZL,phase_velocity./low_freq,d4);
   Zdmax = input_impedance(Zo,ZL,phase_velocity./low_freq,Zmax);
   Zdmin = input_impedance(Zo,ZL,phase_velocity./low_freq,Zmin);

   %%%_Z(f)_OF_STUBS_%%%
   
   Z_PO1 = -i.*Zo.*cot(l_stub_open1.*2.*pi.*low_freq./phase_velocity);
   Z_PO2 = -i.*Zo.*cot(l_stub_open2.*2.*pi.*low_freq./phase_velocity);
   Z_PS1 = i.*Zo.*tan(l_stub_short1.*2.*pi.*low_freq./phase_velocity);
   Z_PS2 = i.*Zo.*tan(l_stub_short2.*2.*pi.*low_freq./phase_velocity);
   Z_SO1 = -i.*Zo.*cot(l_series_open1.*2.*pi.*low_freq./phase_velocity);
   Z_SO2 = -i.*Zo.*cot(l_series_open2.*2.*pi.*low_freq./phase_velocity);
   Z_SS1 = i.*Zo.*tan(l_series_short1.*2.*pi.*low_freq./phase_velocity);
   Z_SS2 = i.*Zo.*tan(l_series_short2.*2.*pi.*low_freq./phase_velocity);
   
   %%%_Z(f)_OF_TUNERS_%%%
   
   ZQmax = input_impedance(ZoPrimeMax,Zdmax,phase_velocity./low_freq,wavelength/4); %QuarterWave at Vmax 
   ZQmin = input_impedance(ZoPrimeMin,Zdmin,phase_velocity./low_freq,wavelength/4); %Quarterwave at Vmin
   ZPO1 = 1./(1./Z_PO1 + 1./Zd1);
   ZPO2 = 1./(1./Z_PO2 + 1./Zd2);
   ZPS1 = 1./(1./Z_PS1 + 1./Zd1);
   ZPS2 = 1./(1./Z_PS2 + 1./Zd2);
   ZSO1 = Z_SO1 + Zd3;
   ZSO2 = Z_SO2 + Zd4;
   ZSS1 = Z_SS1 + Zd3;
   ZSS2 = Z_SS2 + Zd4;
   
   %%%_REFLECTION_COEFFICIENT_OF_TUNERS_%%%
   
   RQmax = reflection(Zo,ZQmax); %QuarterWave at Vmax 
   RQmin = reflection(Zo,ZQmin); %Quarterwave at Vmin
   RPO1 = reflection(Zo,ZPO1);
   RPO2 = reflection(Zo,ZPO2);
   RPS1 = reflection(Zo,ZPS1);
   RPS2 = reflection(Zo,ZPS2);
   RSO1 = reflection(Zo,ZSO1);
   RSO2 = reflection(Zo,ZSO2);
   RSS1 = reflection(Zo,ZSS1);
   RSS2 = reflection(Zo,ZSS2);

   %%%_PLOT_REFLECTION_COEFFICIENT_%%%
   
   plot(low_freq, abs(RPO1));
   hold on
   plot(low_freq, abs(RPO2));
   plot(low_freq, abs(RPS1));
   plot(low_freq, abs(RPS2));
   plot(low_freq, abs(RSO1));
   plot(low_freq, abs(RSO2),'--');
   plot(low_freq, abs(RSS1),'--');
   plot(low_freq, abs(RSS2),'--');
   plot(low_freq, abs(RQmax),'--');
   plot(low_freq, abs(RQmin),'--');
   legend("Open Shunt Solution One", "Open Shunt Solution Two", "Shorted Shunt Solution One", "Shorted Shunt Solution Two", ...
          "Open Series Stub Solution One", "Open Series Stub Solution Two", "Shorted Series Stub Solution One", ...
          "Shorted Series Stub Solution Two", "Quarter Wave Vmax", "Quarter Wave Vmin");
   
   hold off
   title("Reflection Coefficient vs Frequency");
   xlabel("Frequency (Hz)");
   ylabel("Magnitude of Reflection Coefficient");
   
   %%%_FIND_AND_FREQUENCIES_AT_REF=0.2_%%%
   k1 = find(abs(abs(RPO1)-0.2)<0.08);
   k2 = find(abs(abs(RPO2)-0.2)<0.08);
   k3 = find(abs(abs(RPS1)-0.2)<0.08);
   k4 = find(abs(abs(RPS2)-0.2)<0.08);
   k5 = find(abs(abs(RSO1)-0.2)<0.08);
   k6 = find(abs(abs(RSO2)-0.2)<0.08);
   k7 = find(abs(abs(RSS1)-0.2)<0.08);
   k8 = find(abs(abs(RSS2)-0.2)<0.08);
   k9 = find(abs(abs(RQmax)-0.2)<0.08);
   k10 = find(abs(abs(RQmin)-0.2)<0.08);
   
   if (isempty(k1))
       k1 = 1:1:length(low_freq);
   end
   if (isempty(k2))
       k2 = 1:1:length(low_freq);
   end
   if (isempty(k3))
       k3 = 1:1:length(low_freq);
   end
   if (isempty(k4))
       k4 = 1:1:length(low_freq);
   end
   if (isempty(k5))
       k5 = 1:1:length(low_freq);
   end
   if (isempty(k6))
       k6 = 1:1:length(low_freq);
   end
   if (isempty(k7))
       k7 = 1:1:length(low_freq);
   end
   if (isempty(k8))
       k8 = 1:1:length(low_freq);
   end
   if (isempty(k9))
       k9 = 1:1:length(low_freq);
   end
   if (isempty(k10))
       k10 = 1:1:length(low_freq);
   end
   
   band1 = low_freq(k1(length(k1))) - low_freq(k1(1));
   disp("Bandwidth of open shunt solution 1 is: " + band1/(1e9) + " GHz");
   band2 = low_freq(k2(length(k2))) - low_freq(k2(1));
   disp("Bandwidth of open shunt solution 2 is: " + band2/(1e9) + " GHz");
   band3 = low_freq(k3(length(k3))) - low_freq(k3(1));
   disp("Bandwidth of shorted shunt solution 1 is: " + band3/(1e9) + " GHz");
   band4 = low_freq(k4(length(k4))) - low_freq(k4(1));
   disp("Bandwidth of shorted shunt solution 2 is: " + band4/(1e9) + " GHz");
   band5 = low_freq(k5(length(k5))) - low_freq(k5(1));
   disp("Bandwidth of open series stub solution 1 is: " + band5/(1e9) + " GHz");
   band6 = low_freq(k6(length(k6))) - low_freq(k6(1));
   disp("Bandwidth of open series stub solution 2 is: " + band6/(1e9) + " GHz");
   band7 = low_freq(k7(length(k7))) - low_freq(k7(1));
   disp("Bandwidth of shorted series stub solution 1 is: " + band7/(1e9) + " GHz");
   band8 = low_freq(k8(length(k8))) - low_freq(k8(1));
   disp("Bandwidth of shorted series stub solution 2 is: " + band8/(1e9) + " GHz");
   band9 = low_freq(k9(length(k9))) - low_freq(k9(1));
   disp("Bandwidth of Quarter Wave Vmax is: " + band9/(1e9) + " GHz");
   band10 = low_freq(k10(length(k10))) - low_freq(k10(1));
   disp("Bandwidth of Quarter Wave Vmin is: " + band10/(1e9) + " GHz");
   
   
%%%_SMITH_CHART_%%%

   %%%_SETUP_CHART_%%%%

   figure();
   
   th = 0:pi/50:2*pi;
   xunit = cos(th);
   yunit = sin(th);
   plot(xunit,yunit,'k','HandleVisibility','off');
   hold on
   xunit = 0.8.*cos(th) + 0.2;
   yunit = 0.8.*sin(th);
   plot(xunit,yunit,'k','HandleVisibility','off');
   xunit = 0.5.*cos(th) + 0.5;
   yunit = 0.5.*sin(th);
   plot(xunit,yunit,'k','HandleVisibility','off');
   xunit = -1:0.1:1;
   yunit = zeros(1,21);
   plot(xunit,yunit,'k','HandleVisibility','off');
   plot(yunit,xunit,'k','HandleVisibility','off');
   xunit = cos(th) + 1;
   yunit = sin(th) + 1;
   plot(xunit,yunit,'k','HandleVisibility','off');
   plot(xunit,-yunit,'k','HandleVisibility','off');
   xunit = 2.*cos(th) + 1;
   yunit = 2.*sin(th) + 2;
   plot(xunit,yunit,'k','HandleVisibility','off');
   plot(xunit,-yunit,'k','HandleVisibility','off');
   xunit = 5.*cos(th) + 1;
   yunit = 5.*sin(th) + 5;
   plot(xunit,yunit,'k','HandleVisibility','off');
   plot(xunit,-yunit,'k','HandleVisibility','off');
   xlim([-1 1]);
   ylim([-1 1]);   
   pbaspect([1 1 1]);

   %%%_PLOT_DATA_%%%
   plot(RPO1,'LineWidth',3);
   plot(RPO2,'LineWidth',3);
   plot(RPS1,'LineWidth',3);
   plot(RPS2,'LineWidth',3);
   plot(RSO1,'LineWidth',3);
   plot(RSO2,'LineWidth',2);
   plot(RSS1,'LineWidth',2);
   plot(RSS2,'LineWidth',2);
   plot(RQmax,'LineWidth',2);
   plot(RQmin,'LineWidth',2);
   legend("Open Shunt Solution One", "Open Shunt Solution Two", "Shorted Shunt Solution One", "Shorted Shunt Solution Two", ...
          "Open Series Stub Solution One", "Open Series Stub Solution Two", "Shorted Series Stub Solution One", ...
          "Shorted Series Stub Solution Two", "Quarter Wave Vmax", "Quarter Wave Vmin");
      
   plot(real(RPO1(1)),imag(RPO1(1)),'x','HandleVisibility','off');
   plot(real(RPO2(1)),imag(RPO2(1)),'x','HandleVisibility','off');
   plot(real(RPS1(1)),imag(RPS1(1)),'x','HandleVisibility','off');
   plot(real(RPS2(1)),imag(RPS2(1)),'x','HandleVisibility','off');
   plot(real(RSO1(1)),imag(RSO1(1)),'x','HandleVisibility','off');
   plot(real(RSO2(1)),imag(RSO2(1)),'x','HandleVisibility','off');
   plot(real(RSS1(1)),imag(RSS1(1)),'x','HandleVisibility','off');
   plot(real(RSS2(1)),imag(RSS2(1)),'x','HandleVisibility','off');
   plot(real(RQmax(1)),imag(RQmax(1)),'x','HandleVisibility','off');
   plot(real(RQmin(1)),imag(RQmin(1)),'x','HandleVisibility','off');
   
   
   plot(real(RPO1(length(low_freq))),imag(RPO1(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RPO2(length(low_freq))),imag(RPO2(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RPS1(length(low_freq))),imag(RPS1(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RPS2(length(low_freq))),imag(RPS2(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RSO1(length(low_freq))),imag(RSO1(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RSO2(length(low_freq))),imag(RSO2(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RSS1(length(low_freq))),imag(RSS1(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RSS2(length(low_freq))),imag(RSS2(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RQmax(length(low_freq))),imag(RQmax(length(low_freq))),'o','HandleVisibility','off');
   plot(real(RQmin(length(low_freq))),imag(RQmin(length(low_freq))),'o','HandleVisibility','off');
   
   text(0.5,-0.9,'x - start freq, o - stop freq');
   
   hold off
   title("Reflection Coefficient as a function of frequency on Smith Chart");
   xlabel("Real Axis");
   ylabel("Imaginary Axis");
   
%%%_OUTPUT_SOLUTION_INFORMATION_%%%

   %%%_OPEN_SHUNT_SOL_1_%%%
   
   disp("First Solution for Open Shunt Stub");
   disp("    d = " + d1);
   disp("    l = " + l_stub_open1);
   
   %%%_OPEN_SHUNT_SOL_2_%%%
   
   disp("Second Solution for Open Shunt Stub");
   disp("    d = " + d2);
   disp("    l = " + l_stub_open2);
   
   %%%_SHORTED_SHUNT_SOL_1_%%%
   
   disp("First Solution for Shorted Shunt Stub");
   disp("    d = " + d1);
   disp("    l = " + l_stub_short1);
   
   %%%_SHORTED_SHUNT_SOL_2_%%%
   
   disp("Second Solution for Shorted Shunt Stub");
   disp("    d = " + d2);
   disp("    l = " + l_stub_short2);
   
   %%%_OPEN_SERIES_STUB_SOL_1_%%%
   
   disp("First Solution for Open Series Stub");
   disp("    d = " + d3);
   disp("    l = " + l_series_open1);
   
   %%%_OPEN_SERIES_STUB_SOL_2_%%%
   
   disp("Second Solution for Open Series Stub");
   disp("    d = " + d4);
   disp("    l = " + l_series_open2);
   
   %%%_SHORTED_SERIES_STUB_SOL_1_%%%
   
   disp("First Solution for Shorted Series Stub");
   disp("    d = " + d3);
   disp("    l = " + l_series_short1);
   
   %%%_SHORTED_SERIES_STUB_SOL_2_%%%
   
   disp("Second Solution for Shorted Series Stub");
   disp("    d = " + d4);
   disp("    l = " + l_series_short2);
   
   %%%_QUARTER_WAVE_VMAX_SOL_%%%
   
   disp("Quarter Wave Solution at d=Vmax");
   disp("    d = " + Zmax);
   disp("    Zo' = " + ZoPrimeMax);
   
   %%%_QUARTER_WAVE_VMIN_SOL_%%%

   disp("Quarter Wave Solution at d=Vmin");
   disp("    d = " + Zmin);
   disp("    Zo' = " + ZoPrimeMin);
