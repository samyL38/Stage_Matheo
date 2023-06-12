function [new_phase_shift] = mk_same_congruence_phase_shift(phase_shift)
% rearrange phase_shift 


[nb_trials, nb_probes]= size(phase_shift);
new_phase_shift= phase_shift;
for i= 2: nb_trials
    
    ind= abs(new_phase_shift(i,:)- new_phase_shift(i-1,:)) > pi;
    sgn_m= new_phase_shift(i,:)- new_phase_shift(i-1,:) >0;
    sgn_m= sgn_m*(-1);
    sgn_p= new_phase_shift(i,:)- new_phase_shift(i-1,:) <0;
    sgn_modulo= sgn_p + sgn_m;
    new_phase_shift(i,ind)= new_phase_shift(i,ind) + 2*pi*sgn_modulo(ind);
    i;
end