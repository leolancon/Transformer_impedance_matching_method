# Transformer_impedance_matching_method

Script to calculate ideal integrated transformer specifications for millimeter-wave impedance matching

The repository is composed of two MATLAB files:
    1. Transformer_Main.m : This is the main script. Run it to solve impedance matching problems.
    2. TransformerClass.m : This file creates a class used in the main script. Its description is detailed below.
    
For a given impedance matching problem, defined by:
    a source impedance Zs,
    a load impedance Zl,
    an operating frequency f,
    a given coupling coefficient k of the desired integrated transformer,
the script can calculate the self-inductance values of integrated transformers that matches the impedances.


The class TransformerClass contains:
    1. The input of the given matching problem.
    2. The self-inductances of the asymmetrical transformers that can solve the problem.
    3. The equalized sel-inductances of the symmetrical transformers and the discrete components that can solve the problem, following the self-equalizing method.
    4. The required capacitance to add in series with the load or the source if the is no initial solution.
    
    
This class is the base block of the impedance matching method. 
It can be used recursively to create an algorithm that solves impedance matching automatically. This recursive aspect can notably be used to consider the quality factor of the transformer, to limit the value of the self inductances, and to fine-tune the transformers.
