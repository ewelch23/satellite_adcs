# satellite_adcs

This was a project I completed for class. It verifies that the given ADCS requirements can be met with a baseline
design for a satellite.

The requirements are:
### Maneuver 1: Slew 30◦ cross-track.
(a) Perform the maneuver in less than 30 minutes.
(b) The steady state error must be < 0.1◦ .
### Hold attitude constant, i.e. keep the body frame aligned with the inertial frame.
(a) Pointing accuracy must be > 0.15◦
### From an initial state ω = (0.05, 0.05, 0.05) rad/sec, 
bring the spacecraft to ω= (0, 0, 0), ψ = ϕ = θ = 0.
(a) The maneuver must be performed in less than 12 hours.
(b) The steady state error must be < 0.1◦


Adding torque control with a PD controller still does not produce satisfactory answers. 
