This first study is led as an introduction to the fundamental concepts in Finite Element Method
Analisys. For this porpouse it has been taken the example of a truss of beams joined together
by hinges which allow movement only on the plane of the truss and constrained in
all degrees of freedom. The truss is under load imposed to one of the hinges. 

To solve this system and find out the stress to which every beam is subjected it has been 
used a generalized approach to free body diagram: 
- each beam is identified by a number
- each beam extremes (the nodes) are identified by a unic number
- each node degree of freedom is identified by a number

This analitic procedure is the foundation of the FEA method as it allows to solve the 
system using an efficient matrix approach: force and displacement on each node are 
considered as split on two cartesian axis and linked by a linear relationship, whose
coefficient is known as stiffness. Relationships for the mechanical behaviour of the
whole system are also linear and of simple solution. 

The code presented here is an implementation of these powerful concepts to the simple
case proposed.
