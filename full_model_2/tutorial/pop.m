classdef pop < handle
    
    properties
        neurons
        avg_v
        avg_u
        activity
    end
    
    methods
        function obj = pop(E_l, plt)
            neurons = [];
            neurons.append(burster(E_l, plt));
        end
        function r = update(obj,dt)
            for n = obj.neurons
                n.update(dt);
            end
            r = 0;
        end
        function r = draw(obj)
            for n = obj.neurons
                n.draw();
            end
            r = 0;
        end
    end
    
end

