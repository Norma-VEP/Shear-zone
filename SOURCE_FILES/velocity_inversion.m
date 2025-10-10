    

% INVERSION
BC_left     =  BC_left_init;
BC_right    =  BC_right_init;
BC_top      =  BC_top_init;
BC_bottom   =  BC_bottom_init;
if time/SecYear>inversion(1) && time/SecYear<inversion(2)
    BC_left     = BC_left   * (1-(time/SecYear-inversion(1))/(inversion(2)-inversion(1)));
    BC_right    = BC_right  * (1-(time/SecYear-inversion(1))/(inversion(2)-inversion(1)));
    BC_top      = BC_top    * (1-(time/SecYear-inversion(1))/(inversion(2)-inversion(1)));
    BC_bottom   = BC_bottom * (1-(time/SecYear-inversion(1))/(inversion(2)-inversion(1)));
elseif time/SecYear>=inversion(2) && time/SecYear<inversion(3)
    BC_left     =  0;
    BC_right    =  0;
    BC_top      =  0;
    BC_bottom   =  0;
elseif time/SecYear>=inversion(3) && time/SecYear<inversion(4)
    BC_left     = -BC_left   * ((time/SecYear-inversion(3))/(inversion(4)-inversion(3)));
    BC_right    = -BC_right  * ((time/SecYear-inversion(3))/(inversion(4)-inversion(3)));
    BC_top      = -BC_top    * ((time/SecYear-inversion(3))/(inversion(4)-inversion(3)));
    BC_bottom   = -BC_bottom * ((time/SecYear-inversion(3))/(inversion(4)-inversion(3)));
elseif time/SecYear>=inversion(4)
    BC_left     = -BC_left;
    BC_right    = -BC_right;
    BC_top      = -BC_top;
    BC_bottom   = -BC_bottom;
end
if timestep<=1
%     BC_left     =  0;
%     BC_right    =  0;
    BC_top      =  0;
    BC_bottom   =  0;
%     bottom_mix(2,:)  =  0;
end
if timestep>1 && timestep<101
%     BC_left     =  (timestep-21)/10*BC_left;
%     BC_right    =  (timestep-21)/10*BC_right;
    BC_top      =  (timestep-1)/100*BC_top;
    BC_bottom   =  (timestep-1)/100*BC_bottom;
%     bottom_mix(2,:)  =  (timestep-11)/10*bottom_mix(2,:);
end