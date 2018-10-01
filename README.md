# RFR-lapsim
FSAE car lapsim based on 
https://static1.squarespace.com/static/57e8888fc534a547699d733d/t/59705bcf893fc0cc2911099a/1500535765921/WR-217e+Architecture+Design.pdf 
and http://www.jameshakewill.com/Lap_Time_Simulation.pdf

Allows importing a track as a DXF file composed only of lines and arcs, typically exported from a SolidWorks sketch. Make sure that all segments are tangent to their neighbors. 

Simulation takes into account aero downforce and drag, tire friction, and the engine torque curve. Maximum velocities are calculated based on engine output, tire capacbilities, and track characteristics at each point. 

Braking calculations are done by finding the maximum entry velocity for a turn and then back calculating speed at negative relative times based on maximum braking capability. 

Weight transfer is not yet implemented.
