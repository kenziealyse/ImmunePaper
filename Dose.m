function [U] = Dose(t, Dose_amount, Dose_Time)
% Dose_Time./365

if Dose_Time<t && t<Dose_Time+1
    U = Dose_amount;
else
    U = 0;
end

end

