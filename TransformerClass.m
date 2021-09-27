
classdef TransformerClass
    %Theoretical Transformer Derivation
   
    properties
    k                   % Coupling Factor defined by user (unitless)
    freq                % Operating frequency defined by user (GHz)
    Rs                  % Real part of the source impedance defined by user (Ohm)
    Xs                  % Imaginary part of the source impedance (Ohm)
    Xs_user             % Imaginary part of the source impedance defined by user (Ohm)
    Rl                  % Real part of the load impedance defined by user(Ohm)
    Xl                  % Imaginary part of the load impedance (Ohm)
    Xl_user             % Imaginary part of the load impedance defined by user(Ohm)
    omega               % derived from freq
    
    Qs                  % Source Quality factor derived from (Rs,Xs)
    Ql                  % Load Quality factor derived from (Rs,Xs)
    alpha               % normalized parameter used in polynomial root derivations
    Gs                  % Real part of the source admittance derived from Rs & Qs
    Gl                  % Real part of the load admittance derived from Rl & Ql
    
    delta               % Polynomial Determinant
    Zrootp              % Polynomial Positive Root
    Zrootn              % Polynomial Negative Root
    
    Check_sol_rootp_basic      % Variable is 1 if there is a solution for the negative root
    L1_rootp                   % Primary self-inductance of the transformer derived from the positive root
    L2_rootp                   % Secondary self-inductance of the transformer derived from the positive root
    N_rootp                    % Transformer turn ratio derived from the positive root
    Gv_rootp                   % Transformator voltage gain Gv derived from the positive root
    Check_sol_rootn_basic      % Variable is 1 if there is a solution for the negative root
    L1_rootn                   % Primary self-inductance of the transformer derived from the negative root
    L2_rootn                   % Secondary self-inductance of the transformer derived from the negative root
    N_rootn                    % Transformer turn ratio derived from the negative root
    Gv_rootn                   % Transformator voltage gain Gv derived from the negative root
    
    Check_Sol_source_series    % Variable is 1 if there is a solution for Qs source series
    Ltransfo_source_series     % Self-inductance value of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Ltransfo_source_series_1     % Self-inductance value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Ltransfo_source_series_2     % Self-inductance value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Ladj_source_series         % Adjusting series inductor value of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Ladj_source_series_1         % Adjusting series inductor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Ladj_source_series_2         % Adjusting series inductor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Cadj_source_series         % Adjusting series capacitor value of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Cadj_source_series_1         % Adjusting series capacitor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    Cadj_source_series_2         % Adjusting series capacitor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source series adjust
    
    Check_Sol_source_shunt    % Variable is 1 if there is a solution for Qs source shunt
    Ltransfo_source_shunt     % Self-inductance value of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Ltransfo_source_shunt_1     % Self-inductance value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Ltransfo_source_shunt_2     % Self-inductance value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Ladj_source_shunt         % Adjusting series inductor value of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Ladj_source_shunt_1         % Adjusting series inductor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Ladj_source_shunt_2         % Adjusting series inductor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Cadj_source_shunt         % Adjusting series capacitor value of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Cadj_source_shunt_1         % Adjusting series capacitor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    Cadj_source_shunt_2         % Adjusting series capacitor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for source shunt adjust
    
    Check_Sol_load_series    % Variable is 1 if there is a solution for Qs load series
    Ltransfo_load_series     % Self-inductance value of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Ltransfo_load_series_1     % Self-inductance value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Ltransfo_load_series_2     % Self-inductance value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Ladj_load_series         % Adjusting series inductor value of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Ladj_load_series_1         % Adjusting series inductor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Ladj_load_series_2         % Adjusting series inductor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Cadj_load_series         % Adjusting series capacitor value of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Cadj_load_series_1         % Adjusting series capacitor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    Cadj_load_series_2         % Adjusting series capacitor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load series adjust
    
    Check_Sol_load_shunt    % Variable is 1 if there is a solution for Qs load shunt
    Ltransfo_load_shunt     % Self-inductance value of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Ltransfo_load_shunt_1     % Self-inductance value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Ltransfo_load_shunt_2     % Self-inductance value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Ladj_load_shunt         % Adjusting series inductor value of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Ladj_load_shunt_1         % Adjusting series inductor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Ladj_load_shunt_2         % Adjusting series inductor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Cadj_load_shunt         % Adjusting series capacitor value of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Cadj_load_shunt_1         % Adjusting series capacitor value #1 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    Cadj_load_shunt_2         % Adjusting series capacitor value #2 of the symmetrical transformer (L1=L2=Ltransfo) for load shunt adjust
    
    Cnew_max_s                 % additionnal series capacitor on source port if no solution with original user data Xs
    Cnew_max_l                 % additionnal series capacitor on load port if no solution with original user data Xl
    
    end
    
    methods
        function obj = TransformerClass(k,freq,Rs,Xs,Rl,Xl)
        obj.k = k;
        obj.freq = freq*1E9;
        obj.Rs = Rs;
        obj.Xs = Xs;
        obj.Xs_user = Xs;
        obj.Rl = Rl;
        obj.Xl = Xl;
        obj.Xl_user = Xl;
        obj.omega = 2*pi*obj.freq;
        obj.alpha = (1-obj.k^2)/(obj.k^2);
        obj.Qs = -obj.Xs/obj.Rs;
        obj.Ql = -obj.Xl/obj.Rl;
        obj.Gs= 1/(obj.Rs*(1+(obj.Qs)^2));
        obj.Gl= 1/(obj.Rl*(1+(obj.Ql)^2));
        obj.delta=Determinant(obj);
        if obj.delta > 0
            obj=Calcul_L1_L2(obj);
        end
                                            %Interesting fact to note :
                                            %There can be a solution for
                                            %L1=L2 even if the original
                                            %delta < 0.
        obj=Equalize_Source_Series(obj);
        obj=Equalize_Source_Shunt(obj);
        obj=Equalize_Load_Series(obj);
        obj=Equalize_Load_Shunt(obj);
      

        end
        
        

        
       
%-----------------------------------------------------------------------------------------%
 %  Pre-amble : Tools for source or load modification if no solution
 %----------------------------------------------------------------------------------------%       
        
        
        function delta = Determinant (obj)                  % Calculate % Polynomial Determinant
        A = (obj.alpha+1);                                  % Polynomial Coefficient A
        B = -(2*obj.alpha*obj.Qs+obj.Qs+obj.Ql);            % Polynomial Coefficient B
        C = obj.alpha*(obj.Qs^2+1);                         % Polynomial Coefficient C
        delta = B^2-4*A*C;                                 % Polynomial Determinant
        end
        
        function obj = Modification_Zs (obj)  % Calculate the required capacitance to add in series with the source
          delt = Determinant(obj);
          rootn = Calcul_Rootp(obj,obj.Qs,obj.Ql);
          rootp = Calcul_Rootp(obj,obj.Qs,obj.Ql);
             while or(delt <=0, (or(rootn < 0, rootp < 0)))       % Loop to find the minimum impedance shift
                obj.Xs = obj.Xs-1;                                % The additional condition "or(rootn < 0, rootp < 0)" is used to rule out negative self-inductances
                obj.Qs = -obj.Xs/obj.Rs;
                rootp = Calcul_Rootp(obj,obj.Qs,obj.Ql);
                rootn = Calcul_Rootn(obj,obj.Qs,obj.Ql);
                delt = Determinant(obj);
             end
          obj.delta = delt;
          obj.Cnew_max_s = -1/((obj.Xs-obj.Xs_user)*obj.omega);  %Value of the maximum capacitance to add in series to make the matching feasible
          obj.Gs= 1/(obj.Rs*(1+(obj.Qs)^2));
        end
        
        function obj = Modification_Zl (obj)    % Calculate the required capacitance to add in series with the load
          delt = Determinant(obj);
          rootn = Calcul_Rootp(obj,obj.Qs,obj.Ql);
          rootp = Calcul_Rootp(obj,obj.Qs,obj.Ql);
             while or(delt <=0, (or(rootn < 0, rootp < 0)))                     % Loop to find the minimum impedance shift
                obj.Xl = obj.Xl-1;                                              % The additional condition "or(rootn < 0, rootp < 0)" is used to rule out negative self-inductances
                obj.Ql = -obj.Xl/obj.Rl;
                rootp = Calcul_Rootp(obj,obj.Qs,obj.Ql);
                rootn = Calcul_Rootn(obj,obj.Qs,obj.Ql);
                delt = Determinant(obj);
             end
          obj.delta = delt;
          obj.Cnew_max_l = -1/((obj.Xl-obj.Xl_user)*obj.omega);     %Value of the maximum capacitance to add in series to make the matching feasible
          obj.Gl= 1/(obj.Rl*(1+(obj.Ql)^2));                        %This value is given for reference only. A new object of class Transformer must be created  
        end                                                         %with the corrected impedance (including the added cap.) to take it into account.
        
        
%-----------------------------------------------------------------------------------------%
 % I. Derivation for the general case of a non-symetrical Transformer :  (L1,L2) different
 %----------------------------------------------------------------------------------------%
        
        function obj = Calcul_L1_L2(obj)        % Calculate the theoretical transformer (L1,L2) for the positive root (rootp) and for the negative root (rootn)
            obj.Zrootn = Calcul_Rootn (obj,obj.Qs,obj.Ql);     % Polynomial Negative Root (Zrootn)
            obj.Zrootp = Calcul_Rootp (obj,obj.Qs,obj.Ql);     % Polynomial Positive Root (Zrootp)   
            obj.L1_rootp = obj.Zrootp * obj.Rs/(obj.k^2*obj.omega*obj.alpha);                               %Solution for L1 with positive root
            obj.L2_rootp = obj.k^2*obj.L1_rootp*obj.Rl*(1+obj.Ql^2)/(obj.Rs+obj.Rs*(obj.Qs-obj.Zrootp)^2);  %Solution for L2 with positive root
            obj.N_rootp = obj.k*(obj.L1_rootp/obj.L2_rootp)^0.5;
            obj.Gv_rootp = 1/abs((obj.N_rootp*(1+obj.alpha+1i*obj.alpha*obj.omega*obj.L2_rootp/(obj.Rl+1i*obj.Xl_user))));
            obj.L1_rootn = obj.Zrootn*obj.Rs/(obj.k^2*obj.omega*obj.alpha);                                 %Solution for L1 with negative root
            obj.L2_rootn = obj.k^2*obj.L1_rootn*obj.Rl*(1+obj.Ql^2)/(obj.Rs+obj.Rs*(obj.Qs-obj.Zrootn)^2);  %Solution for L2 with positive root
            obj.N_rootn = obj.k*(obj.L1_rootn/obj.L2_rootn)^0.5;
            obj.Gv_rootn = 1/abs((obj.N_rootn*(1+obj.alpha+1i*obj.alpha*obj.omega*obj.L2_rootn/(obj.Rl+1i*obj.Xl_user))));
        
            if (obj.L1_rootp > 1E-12 && obj.L2_rootp > 1E-12)  %Check if the found self-inductance values make sense.
                obj.Check_sol_rootp_basic = 1;
            else
                obj.Check_sol_rootp_basic = 0;
                obj.L1_rootp = 0;                              
                obj.L2_rootp = 0;  
                obj.N_rootp = 0;
            end
        
            if (obj.L1_rootn > 1E-12 && obj.L2_rootn > 1E-12) %Check if the found self-inductance values make sense.
                obj.Check_sol_rootn_basic = 1;
            else
                obj.Check_sol_rootn_basic = 0;
                obj.L1_rootn = 0;                              
                obj.L2_rootn = 0;  
                obj.N_rootn = 0;
            end
            
            
        end
            
        function rootn = Calcul_Rootn (obj,Qs,Ql)    % Calculate the negative root (rootn) for the 2nd polynomial derived from transformer-based (power) matching conditions
        A = (obj.alpha+1);            % Polynomial Coefficient A
        B = -(2*obj.alpha*Qs+Qs+Ql);  % Polynomial Coefficient B
        C = obj.alpha*(Qs^2+1);       % Polynomial Coefficient C
        delta = B^2-4*A*C;            % Polynomial Determinant
        rootn = (-B-delta^.5)/(2*A);  % Polynomial Negative Root
        end
        
        function rootp = Calcul_Rootp (obj,Qs,Ql)    % Calculate the positive root (rootp) for the 2nd polynomial derived from transformer-based (power) matching conditions
        A = (obj.alpha+1);            % Polynomial Coefficient A
        B = -(2*obj.alpha*Qs+Qs+Ql);  % Polynomial Coefficient B
        C = obj.alpha*(Qs^2+1);       % Polynomial Coefficient C
        delta = B^2-4*A*C;            % Polynomial Determinant
        rootp = (-B+delta^.5)/(2*A);  % Polynomial Negative Root
        end
        
 %----------------------------------------------------------------------------------------%
 % II. Derivation for Symmetrical Transformer : L1=L2
 %----------------------------------------------------------------------------------------%

        %------------- Preliminary Functions ------------%
 
        function [is_inductor,LC] = Calcul_Element_Adjust (obj,Delta_Q,is_series,is_source)
            
            if (is_series)
                if (is_source)
                    R=obj.Rs;
                else
                    R=obj.Rl;
                end
                if Delta_Q < 0    % Sign of Delta_Qs gives the type of compensating component (C1 or L1)
                    is_inductor = true ;
                    LC = R * abs(Delta_Q) / (obj.omega);
                else
                    is_inductor = false;
                    LC = 1 / (abs(Delta_Q)* R * obj.omega);
                end
            else
                if (is_source)
                    G=obj.Gs;
                else
                    G=obj.Gl;
                end
                if Delta_Q < 0    % Sign of Delta_Qs gives the type of compensating device (C1 or L1)
                   is_inductor = true ;
                   LC = 1 / (abs(Delta_Q)* G * obj.omega);
                else
                  is_inductor = false;
                  LC = (abs(Delta_Q)* G) / (obj.omega);
               end
            end   
        end
            
 
 
        %------------- II.A/ Source Adjust (Qs) with Series Device ------------%
        function obj = Equalize_Source_Series (obj)                          %Workout the required compensation device on port 1 (C1 or L1) to enable L1=L2=L_rootp_port1 for the positive root
            obj.Check_Sol_source_series = false;
            is_series = true;
            is_source = true;
            Element=[];
            Leq=[];
            
            if (obj.k^2*(obj.Rl*(1+obj.Ql^2)/obj.Rs)-1) >= 0 
                A = sqrt(obj.k^2*(obj.Rl*(1+obj.Ql^2)/obj.Rs)-1);
                Qs_pos_new=(A*(A*(obj.alpha+1)+obj.Ql)+obj.alpha)/(obj.Ql+A);      %Derive a Qs solution, named Qs_pos_new, associated to one root (n or p we don't know so far)
                Qs_neg_new=(A*(A*(obj.alpha+1)-obj.Ql)+obj.alpha)/(obj.Ql-A);      %Derive an other Qs,named Qs_neg_new, associated to the other root (n or p we don't know so far)
                
                Qs_list = [Qs_pos_new,Qs_neg_new]; %List of the new Qs that allows L1=L2
                
                
                for idx = 1:length(Qs_list)
                    q = Qs_list(idx);
                    if (abs(obj.Rs/(1-obj.k^2)-obj.Rl/obj.alpha*(1+obj.Ql^2)/(1+(q-Calcul_Rootp(obj,q,obj.Ql))^2))<1E-3) % this test allows for recovering which Qs (Qs_pos_new or Qs_neg_new) is actually linked to the positive root (rootp)
                        Delta_Qs = q - obj.Qs;
                        [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Qs,is_series,is_source);
                        Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                        Zrootp_new = Calcul_Rootp(obj,q,obj.Ql);                                           % The new polynomial root is calculated
                        Leq =  [Leq;Zrootp_new*obj.Rs/(obj.k^2*obj.omega*obj.alpha)];                           % vector of equivalent transformer values with negative root                                                                    % vector of equivalent transformer values with positive root
                    else
                        if (abs(obj.Rs/(1-obj.k^2)-obj.Rl/obj.alpha*(1+obj.Ql^2)/(1+(q-Calcul_Rootn(obj,q,obj.Ql))^2))<1E-3) % this test allows for recovering which Qs (Qs_pos_new or Qs_neg_new) is actually linked to the negative root (rootn)
                            Delta_Qs = q - obj.Qs;
                            [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Qs,is_series,is_source);
                            Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with negative root
                            Zrootn_new = Calcul_Rootn(obj,q,obj.Ql);                                           % The new polynomial root is calculated
                            Leq =  [Leq;Zrootn_new*obj.Rs/(obj.k^2*obj.omega*obj.alpha)];                            % vector of equivalent transformer values with negative root
                        end
                    end
                    
                    if (Leq(idx)>1E-12 && Element(idx,2)>1E-15 )
                        obj.Check_Sol_source_series = true;
                        obj.Ltransfo_source_series(idx) = Leq(idx);
                        if (Element(idx,1))
                            obj.Ladj_source_series(idx) = Element(idx,2);
                            obj.Cadj_source_series(idx) = 0;
                        else
                            obj.Ladj_source_series(idx) = 0;
                            obj.Cadj_source_series(idx) = Element(idx,2);
                        end
                    end
                end
            end
            
            if ((obj.Check_Sol_source_series == true) && (length(obj.Ltransfo_source_series) == 1))
            obj.Check_Sol_source_series = 1;
            obj.Ltransfo_source_series_1 = obj.Ltransfo_source_series(1);
            obj.Ladj_source_series_1 = obj.Ladj_source_series(1);
            obj.Cadj_source_series_1 = obj.Cadj_source_series(1);
            obj.Ltransfo_source_series_2 = 0;
            obj.Ladj_source_series_2 = 0;
            obj.Cadj_source_series_2 = 0;

           else
               if ((obj.Check_Sol_source_series == true) && (length(obj.Ltransfo_source_series) == 2))
                    obj.Check_Sol_source_series = 2; 
                    obj.Ltransfo_source_series_1 = obj.Ltransfo_source_series(1);
                    obj.Ladj_source_series_1 = obj.Ladj_source_series(1);
                    obj.Cadj_source_series_1 = obj.Cadj_source_series(1);
                    obj.Ltransfo_source_series_2 = obj.Ltransfo_source_series(2);
                    obj.Ladj_source_series_2 = obj.Ladj_source_series(2);
                    obj.Cadj_source_series_2 = obj.Cadj_source_series(2);
               end
           end
   
        end

                
        %------------- II.B/ Source Adjust (Qs) with Shunt Device ------------%
        
        function obj = Equalize_Source_Shunt (obj)                    % Workout the required compensation device on port 1 (C1 or L1) to enable L1=L2=L_rootn_port1 for the negative root
            obj.Check_Sol_source_shunt = false;
            is_series = false;
            is_source = true;
            a = 4*(obj.alpha+1)^2/(obj.alpha*obj.Gl)*obj.Gs*(1-obj.k^2)*((obj.alpha+1)^2/(obj.alpha*obj.Gl)*obj.Gs*(1-obj.k^2)-1);
            b = -8*obj.Ql*obj.Gs*(obj.alpha+1)^2/obj.Gl*(1-obj.k^2);
            c = 4*((obj.alpha+1)^2*(1+obj.Ql^2)+(obj.alpha+1)^2/(obj.alpha*obj.Gl)*obj.Gs*(1-obj.k^2)*(2*(obj.alpha+1)^2/(obj.alpha*obj.Gl)*obj.Gs*(1-obj.k^2)-obj.Ql^2-2*obj.alpha-3));
            d = b;
            e = (2*(obj.alpha+1)^2/(obj.alpha*obj.Gl)*obj.Gs*(1-obj.k^2)-obj.Ql^2-2*(obj.alpha+1))^2-obj.Ql^4+4*obj.alpha*obj.Ql^2*(obj.alpha+1);
            
            
            p=[a b c d  e]; % define a vector featuring the polynomial coefficient
            r=roots(p); % derive the roots of the polynomial featuring the coefficients defined in vector p
            Qs_list=[]; % list of new Qs that achieve L1=L2 for the case
            
            Element=[];
            Leq=[];
            
            for idx=1:length(r)
                if r(idx)==real(r(idx))
                    Qs_list=[Qs_list;r(idx)];
                end
            end
            
            for idx=1:length(Qs_list)
                q=Qs_list(idx);
                if (abs(obj.Gs*(1+q^2)*(1-obj.k^2)-obj.alpha*obj.Gl*(1+(q-Calcul_Rootn(obj,q,obj.Ql))^2))<1E-3) % this test allows for recovering which Qs is actually linked to the negative root (rootn)
                    Delta_Qs = q - obj.Qs;
                    [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Qs,is_series,is_source);
                    Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                    Zrootn_new = Calcul_Rootn(obj,q,obj.Ql);                                           % The new polynomial root is calculated
                    Leq =  [Leq; Zrootn_new/((1-obj.k^2)*(1+q^2)*obj.omega*obj.Gs)];                           % vector of equivalent transformer values with negative root                                                                    % vector of equivalent transformer values with positive root
                else
                    if (abs(obj.Gs*(1+q^2)*(1-obj.k^2)-obj.alpha*obj.Gl*(1+(q-Calcul_Rootp(obj,q,obj.Ql))^2))<1E-3) % this test allows for recovering which Qs is actually linked to the positive root (rootp)
                        Delta_Qs = q - obj.Qs;
                        [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Qs,is_series,is_source);
                        Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                        Zrootp_new = Calcul_Rootp(obj,q,obj.Ql);                                           % The new polynomial root is calculated
                        Leq =  [Leq;Zrootp_new/((1-obj.k^2)*(1+q^2)*obj.omega*obj.Gs)];                            % vector of equivalent transformer values with positive root
                    end
                end
                
                if (Leq(idx)>1E-12 && Element(idx,2)>1E-15 )
                    obj.Check_Sol_source_shunt = true;
                    obj.Ltransfo_source_shunt(idx) = Leq(idx);
                    if (Element(idx,1))
                        obj.Ladj_source_shunt(idx) = Element(idx,2);
                        obj.Cadj_source_shunt(idx) = 0;
                    else
                        obj.Ladj_source_shunt(idx) = 0;
                        obj.Cadj_source_shunt(idx) = Element(idx,2);
                    end
                end
            end
            
            if ((obj.Check_Sol_source_shunt == true) && (length(obj.Ltransfo_source_shunt) == 1))
            obj.Check_Sol_source_shunt = 1;
            obj.Ltransfo_source_shunt_1 = obj.Ltransfo_source_shunt(1);
            obj.Ladj_source_shunt_1 = obj.Ladj_source_shunt(1);
            obj.Cadj_source_shunt_1 = obj.Cadj_source_shunt(1);
            obj.Ltransfo_source_shunt_2 = 0;
            obj.Ladj_source_shunt_2 = 0;
            obj.Cadj_source_shunt_2 = 0;
            else
               if ((obj.Check_Sol_source_shunt == true) && (length(obj.Ltransfo_source_shunt) == 2))
                    obj.Check_Sol_source_shunt = 2; 
                    obj.Ltransfo_source_shunt_1 = obj.Ltransfo_source_shunt(1);
                    obj.Ladj_source_shunt_1 = obj.Ladj_source_shunt(1);
                    obj.Cadj_source_shunt_1 = obj.Cadj_source_shunt(1);
                    obj.Ltransfo_source_shunt_2 = obj.Ltransfo_source_shunt(2);
                    obj.Ladj_source_shunt_2 = obj.Ladj_source_shunt(2);
                    obj.Cadj_source_shunt_2 = obj.Cadj_source_shunt(2);
               end
            end
       
        end
        
        %------------- II.C/ Load Adjust (Ql) with Series Device --------------%
        
        function obj = Equalize_Load_Series (obj)                          %Workout the required compensation device on port 2 (C2 or L2) to enable L1=L2=L_rootn_port2 with a series element
            obj.Check_Sol_load_series = false;
            is_series = true;
            is_source = false;
            
            gamma = obj.k^2*obj.Rl/obj.Rs;
            a=gamma*(gamma*(obj.alpha+1)^2-1);
            b=-2*obj.alpha*obj.Qs*gamma;
            c=2*gamma*(gamma-1)*(obj.alpha+1)^2+1+obj.Qs^2-gamma*(1+obj.Qs^2-2*obj.alpha*(obj.alpha+1));
            d=b;
            e=(obj.alpha+1)^2*(gamma*(gamma-2)+1)+(1-gamma)*(obj.Qs^2-2*obj.alpha*(obj.alpha+1))+obj.alpha^2;
            
            p=[a b c d e]; % define a vector featuring the polynomial coefficient
            r=roots(p); % derive the roots of the polynomial featuring the coefficients defined in vector p
            Ql_list=[]; % list of QL new that achieve L1=L2 for the case
            
            Element=[];
            Leq=[];
            
            for idx=1:length(r)
                if r(idx)==real(r(idx))
                    Ql_list=[Ql_list;r(idx)];
                end
            end
            
            for idx=1:length(Ql_list)
                q=Ql_list(idx);
                if (abs(obj.Rs/(1-obj.k^2)-obj.Rl/obj.alpha*(1+q^2)/(1+(obj.Qs-Calcul_Rootn(obj,obj.Qs,q))^2))<1E-3) % this test allows for recovering which Qs (Qs_pos_new or Qs_neg_new) is actually linked to the negative root (rootn)
                    Delta_Ql = q - obj.Ql;
                    [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Ql,is_series,is_source);
                    Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                    Zrootn_new = Calcul_Rootn(obj,obj.Qs,q);                                           % The new polynomial root is calculated
                    Leq =  [Leq; Zrootn_new*obj.Rs/(obj.k^2*obj.omega*obj.alpha)];                           % vector of equivalent transformer values with negative root                                                                    % vector of equivalent transformer values with positive root
                else
                    if (abs(obj.Rs/(1-obj.k^2)-obj.Rl/obj.alpha*(1+q^2)/(1+(obj.Qs-Calcul_Rootp(obj,obj.Qs,q))^2))<1E-3) % this test allows for recovering which Qs is actually linked to the positive root (rootp)
                        Delta_Ql = q - obj.Ql;
                        [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Ql,is_series,is_source);
                        Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                        Zrootp_new = Calcul_Rootp(obj,obj.Qs,q);                                           % The new polynomial root is calculated
                        Leq =  [Leq; Zrootp_new*obj.Rs/(obj.k^2*obj.omega*obj.alpha)];                            % vector of equivalent transformer values with positive root
                    end
                end
                
                if (Leq(idx)>1E-12 && Element(idx,2)>1E-15 )
                    obj.Check_Sol_load_series = true;
                    obj.Ltransfo_load_series(idx) = Leq(idx);
                    if (Element(idx,1))
                        obj.Ladj_load_series(idx) = Element(idx,2);
                        obj.Cadj_load_series(idx) = 0;
                    else
                        obj.Ladj_load_series(idx) = 0;
                        obj.Cadj_load_series(idx) = Element(idx,2);
                    end
                end
            end
            
            if ((obj.Check_Sol_load_series == true) && (length(obj.Ltransfo_load_series) == 1))
            obj.Check_Sol_load_series = 1;
            obj.Ltransfo_load_series_1 = obj.Ltransfo_load_series(1);
            obj.Ladj_load_series_1 = obj.Ladj_load_series(1);
            obj.Cadj_load_series_1 = obj.Cadj_load_series(1);
            obj.Ltransfo_load_series_2 = 0;
            obj.Ladj_load_series_2 = 0;
            obj.Cadj_load_series_2 = 0;
            else
               if ((obj.Check_Sol_load_series == true) &&  (length(obj.Ltransfo_load_series) == 2))
                    obj.Check_Sol_load_series = 2; 
                    obj.Ltransfo_load_series_1 = obj.Ltransfo_load_series(1);
                    obj.Ladj_load_series_1 = obj.Ladj_load_series(1);
                    obj.Cadj_load_series_1 = obj.Cadj_load_series(1);
                    obj.Ltransfo_load_series_2 = obj.Ltransfo_load_series(2);
                    obj.Ladj_load_series_2 = obj.Ladj_load_series(2);
                    obj.Cadj_load_series_2 = obj.Cadj_load_series(2);
               end
            end
            
        end
  
       %------------- II.D/ Load Adjust (Ql) with Shunt Device ------------% 
        
       function obj = Equalize_Load_Shunt (obj)                          %Workout the required compensation device on port 2 (C2 or L2) to enable L1=L2=L_rootn_port2 with a shunt element
           
           obj.Check_Sol_load_shunt = false;
           is_series = false;
           is_source = false;
           
           C=obj.Gs*(1+obj.Qs^2)*(1-obj.k^2)/(obj.alpha*obj.Gl)-1;
           a=C*(obj.alpha+1)^2-obj.Qs^2*(1+obj.alpha*(obj.alpha+2));
           b=2*obj.alpha*obj.Qs*(obj.alpha+1)*(1+obj.alpha+C*(obj.alpha+1));
           c=(obj.alpha+1)^2*(C*obj.Qs^2-(obj.alpha+C*(obj.alpha+1))^2);
           
           p=[a b c ]; % define a vector featuring the polynomial coefficient
           r=roots(p); % derive the roots of the polynomial featuring the coefficients defined in vector p
           Ql_list=[]; % list of QL new that achieve L1=L2 for the case
           
           Element=[];
           Leq=[];
           
           for idx=1:length(r)
               if r(idx)==real(r(idx))
                   Ql_list=[Ql_list;r(idx)];
               end
           end
           
           
           for idx=1:length(Ql_list)
               q=Ql_list(idx);
               if (abs(obj.Gs*(1+obj.Qs^2)*(1-obj.k^2)-obj.alpha*obj.Gl*(1+(obj.Qs-Calcul_Rootn(obj,obj.Qs,q))^2))<1E-3) % this test allows for recovering which Ql (Qs_pos_new or Qs_neg_new) is actually linked to the negative root (rootn)
                   Delta_Ql = q - obj.Ql;
                   [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Ql,is_series,is_source);
                   Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with negative root
                   Zrootn_new = Calcul_Rootn(obj,obj.Qs,q);                                           % The new polynomial root is calculated
                   Leq =  [Leq; Zrootn_new/((1-obj.k^2)*(1+obj.Qs^2)*obj.omega*obj.Gs)];                           % vector of equivalent transformer values with negative root                                                                    % vector of equivalent transformer values with positive root
               else
                   if (abs(obj.Gs*(1+obj.Qs^2)*(1-obj.k^2)-obj.alpha*obj.Gl*(1+(obj.Qs-Calcul_Rootp(obj,obj.Qs,q))^2))<1E-3) % this test allows for recovering which Qs is actually linked to the positive root (rootp)
                       Delta_Ql = q - obj.Ql;
                       [is_inductor,LC]=Calcul_Element_Adjust(obj,Delta_Ql,is_series,is_source);
                       Element = [Element;[is_inductor,LC]];   % vector of equivalent adjust element with positive root
                       Zrootp_new = Calcul_Rootp(obj,obj.Qs,q);                                           % The new polynomial root is calculated
                       Leq =  [Leq; Zrootp_new/((1-obj.k^2)*(1+obj.Qs^2)*obj.omega*obj.Gs)];                            % vector of equivalent transformer values with positive root
                   end
               end
               
               if (Leq(idx)>1E-12 && Element(idx,2)>1E-15 )  
                   obj.Check_Sol_load_shunt = true;
                   obj.Ltransfo_load_shunt(idx) = Leq(idx);
                   if (Element(idx,1))
                       obj.Ladj_load_shunt(idx) = Element(idx,2);
                       obj.Cadj_load_shunt(idx) = 0;
                   else
                       obj.Ladj_load_shunt(idx) = 0;
                       obj.Cadj_load_shunt(idx) = Element(idx,2);
                   end
               end
           end
           
           if ((obj.Check_Sol_load_shunt == true) && (length(obj.Ltransfo_load_shunt) == 1))
            obj.Check_Sol_load_shunt = 1;
            obj.Ltransfo_load_shunt_1 = obj.Ltransfo_load_shunt(1);
            obj.Ladj_load_shunt_1 = obj.Ladj_load_shunt(1);
            obj.Cadj_load_shunt_1 = obj.Cadj_load_shunt(1);
            obj.Ltransfo_load_shunt_2 = 0;
            obj.Ladj_load_shunt_2 = 0;
            obj.Cadj_load_shunt_2 = 0;
           else
               if ((obj.Check_Sol_load_shunt == true) && (length(obj.Ltransfo_load_shunt) == 2))
                    obj.Check_Sol_load_shunt = 2; 
                    obj.Ltransfo_load_shunt_1 = obj.Ltransfo_load_shunt(1);
                    obj.Ladj_load_shunt_1 = obj.Ladj_load_shunt(1);
                    obj.Cadj_load_shunt_1 = obj.Cadj_load_shunt(1);
                    obj.Ltransfo_load_shunt_2 = obj.Ltransfo_load_shunt(2);
                    obj.Ladj_load_shunt_2 = obj.Ladj_load_shunt(2);
                    obj.Cadj_load_shunt_2 = obj.Cadj_load_shunt(2);
               end
           end
       end
       
    end




        
end