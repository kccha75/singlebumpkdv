function RHS=findRHS(x,gamma)

% Define RHS as in thesis
% Equation Parameters (see notes)
l=0.5;
RHS=zeros(length(x),1);

for i=1:length(x)
    if abs(x(i))<=l
        RHS(i)=(4/(3*l))*cos(pi/2*x(i)./l)^4;
    else
        RHS(i)=0;
    end
end

RHS=RHS*gamma;

end