classdef World
    properties
        L;
        v_max;
        ni;nj;
        dx;dv;
        periodic = true;
    end
    methods
        %         function obj = World()
        %         end
        function obj = setLimits(obj,L,v_max)
            obj.L = L;
            obj.v_max = v_max;
        end
        function obj = setNodes(obj,N,M)
            obj.ni = N;
            obj.nj = 2*M-1;
            obj.dx = obj.L/(obj.ni-1);
            obj.dv = 2*obj.v_max/(obj.nj-1);
        end
        function x = getX(obj,i)
            x = 0+(i-1)*obj.dx;
        end
        function v = getV(obj,j)
            v = -obj.v_max+(j-1)*obj.dv;
        end

        %linear interpolation: a higher order scheme needed!
        function val = interp(obj,f,x,v)
            fi = x/obj.dx+1;
            fj = (v-(-obj.v_max))/obj.dv+1;

            if(obj.periodic)
                if(fi<1)
                    fi = fi+obj.ni-1;
                end
                if(fi>obj.ni)
                    fi = fi-obj.ni+1;
                end

            else,if(fi<1 || fi>=obj.ni)
                    val = 0;
                    return
                end
            end

            if(fj<1||fj>obj.nj)
                val = 0;
                return
            end

            i = fix(fi);
            j = fix(fj);
            di = fi-i;
            dj = fj-j;

            val = (1-di)*(1-dj)*f(i,j);
            if(i<obj.ni)
                val = val+di*(1-dj)*f(i+1,j);
            end
            if(j<obj.nj)
                val = val+(1-di)*dj*f(i,j+1);
            end
            if(i<obj.ni && j<obj.nj)
                val = val+di*dj*f(i+1,j+1);
            end
        end

        % makes values on left and right edge identical on periodic systems
        function f = applyBC(obj,f)
            if(~obj.periodic)
                return;
            end
            for j = 1:obj.nj
                f(1,j) = 0.5*(f(1,j)+f(obj.ni,j));
                f(obj.ni,j) = f(1,j);
            end
        end

%         function data = filter(~,a)
%             if a<abs(1e-20)
%                 data = 0;
%             else
%                 data = a;
%             end
%         end
    end
end
