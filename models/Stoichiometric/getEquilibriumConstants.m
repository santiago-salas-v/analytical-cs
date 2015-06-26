function [ output_args ] = getEquilibriumConstants( input_args )
%GETEQUILIBRIUMCONSTANTS Calculates equilibrium constants at T=298.15K
%   inputs include reaction coefficient matrix. Ex. reaction set 
%   A <<==>> B, B ==>> C is written as [-1,+1,0;+1,-1,0;0,-1,+1]
[~,standard_thermodynamic_properties]  ...
                = cargarCSV('standard_thermodynamic_properties.csv');

end