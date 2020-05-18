%%%_Jonathan_and_Jodlet_%%%
function z = input_impedance(Zo, ZL, wavelength, L)
    B = 2*pi./wavelength;
    z = Zo.*(ZL+i.*Zo.*tan(B.*L))./(Zo + i.*ZL.*tan(B.*L));
end