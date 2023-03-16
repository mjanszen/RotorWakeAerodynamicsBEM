# BEM code - Rotor wake 


## Structure of the code 

1. Read in input data (Airfoils, twist, chord, etc...)
2. Initialize induction
3. Calculate velocity and loads at the blade
4. Calculate new induction
5. test convergence ---> restart from 2 
6. When induction is converged, calculate all following values 


### To dos: 
- copy functions from old codes 
- make a general code setup 
- define functions and interface information between functions   


- function input values 

- Function Cl, CD interpolate from Airfoil polar (In: aoa [rad]; out [Cl, CD]  ) 
  - 
- Function Ct Cn (IN: AOA, Cl, CD; OUT: [Ct, Cn])

- Function for induction factors (In: r,V, ...  )
  - Prandtl tip loss fct (IN: radial position, phi,n_blades ; OUT: tip loss factor)
  - Glauert correction 
- Function for input 

- plot functions 
- stay sexy