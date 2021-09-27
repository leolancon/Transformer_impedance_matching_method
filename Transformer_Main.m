%          Impedance matching method using integrated transformers
%                   Thierry Taris - Léo Lançon
%-----------------------------------------------------------------------

%------- Script description--------

%This script allow to calculate the required characteristics of an
%integrated transformer to match a given set of impedances.

% Read the associated README in the github repository for more details. 

%-----------Script-----------------



% Specifications - Set up by the users in the "input" vector. Can be done
% in the command window or directly in the script.



Spec = [];

Spec = input('Define: \n k freq(GHz) Rs(Ohm) Xs(Ohm) Rl(Ohm) Xl(Ohm) \n');   % Setup inputs as a vector (ex: [0.8 77 100 -200 50 -100], 
                                                                             % to adapt Zs=100-j200 to Zl=50-j100, at 77 GHz, with a coupling coefficient of 0.8
                                                                             
T = TransformerClass(Spec(1),Spec(2),Spec(3),Spec(4),Spec(5),Spec(6)) % Create an object of class TransformerBis (see TransformerBis.m), 
                                                                      % which will contain the existing solutions for the given inputs

                                                                      
                                                                      
% The following part of the script can be used in case the
% matching problem does not have a direct solution.

% If so, uncomment the following code

% delta=T.Determinant();
% if or(delta <=0, (or(T.Zrootn < 0, T.Zrootp < 0)))                                
%     disp('no solution');
%     Port_Choice = input('Which root? \n 1:source or 2:load \n');
%     Select_port = Port_Choice(1);
% 
%     if Select_port == 1
%         T = T.Modification_Zs();
%     else
%         T = T.Modification_Zl();
%     end
%     T
%         
% else
%         
% end    
