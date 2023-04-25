function [Vbus, Iline, Iload] = convert_to_single_vector(Vbus_3Phase, Iline_3Phase, Iload_3Phase, Sys, branch_number)
% This function converts the three-phase voltage, branch current, and load current 
% vectors into a single vector format for each of them.

Vbus=zeros(3*Sys.bus_number,1);
for k=1:1:Sys.bus_number
  for k1= 1:1:3
    Vbus((k-1)*3+k1)=Vbus_3Phase(k1,1,k);
  end
end

% change three phase branch current to a single vector
Iline=zeros(3*branch_number,1);
for k=1:1:branch_number
  for k1= 1:1:3
    Iline((k-1)*3+k1)=Iline_3Phase(k1,1,k);
  end
end

% change three phase load current to a single vector
Iload=zeros(3*length(Sys.LoadIdx),1);
for k=1:1:length(Sys.LoadIdx)
  for k1= 1:1:3
    Iload((k-1)*3+k1)=Iload_3Phase(k1,1,k);
  end
end   
