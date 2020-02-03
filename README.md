# Nuclear growth
A simulation of nuclear growth from pER.

## Usage
Start matlab, change directory to where this code is. Then run
```matlab
startup
```
To add the functions to the path.  

To simulate the nuclear growth with defaults parameters, run :
```matlab
script_integ_nucleus
```

To find the best set of parameters for nuclear growth sequence, run :
```matlab
script_fit_nucleus
```

## Changing parameters
To change parameters, you can modify them in a script calling nuclear_growth_sequence, e.g. :  
```matlab
...
parameters=nucleus_parameters();
parameters.N_sat=parameters.N_sat/2;
[CELL_STAGE,flag]=nuclear_growth_sequence(exp_data,parameters);
```

Or create a new parameter file, e.g. : parameters/altered_parameters.m  
In which case the first line of parameters/altered_parameters.m should be :
```matlab
function [ parameters ] = altered_parameters(  )
```

And then call it in your script (e.g. copied from script_integ_nucleus) as :  
```matlab
...
parameters=altered_parameters();
...
[CELL_STAGE,flag]=nuclear_growth_sequence(exp_data,parameters);
```

## Serge Dmitrieff -- http://biophysics.fr
