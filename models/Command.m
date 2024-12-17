classdef Command < Simulink.IntEnumType
    %COMMAND enumeration
    %   list of command received from serial communication
    % command format:
    % [Mode_Code, Parameter1, Parameter2]
    
    enumeration
        STOP(0)
        SWITCHON(10)
        SHUTDOWN(11)
        ENABLE(12)
        DISABLE(13)
        POSITIONMODE(14)
        EXITPOSITIONMODE(15)
        SPEEDMODE(16)
        EXITSPEEDMODE(17)

        POSITIONCOMMAND(20) % + parameter1 + parameter2
        SPEEDCOMMAND(21) % + parameter1 + parameter2
    end
end

% Status code

Ready to switch on: 50
Switched on: 51
Operation enable: 52
POSITIONMODE: 53
