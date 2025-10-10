alpha    =   wm*dt;
% Rotate stress on marker according to its spin
% Compute amount of rotation from spin rate:
% Espin=1/2(dvy/dx-dvx/dy) i.e. positive for clockwise rotation
% (when x axis is directed rightward and y axis is directed downward)

Txxm = Txxm + dtxxm;
Txym = Txym + dtxym;

txym    =   Txxm.*sin(2.*alpha) + Txym.*cos(2.*alpha);
txxm    =   Txxm.*((cos(alpha)).^2) - Txxm.*((sin(alpha)).^2) - Txym.*sin(2.*alpha);

Txxm = txxm;
Txym = txym;

%==============================================================
%==========================================================
