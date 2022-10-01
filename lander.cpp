// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"

vector3d drag_force(vector3d &v, vector3d &pos) {
    double CdA = DRAG_COEF_LANDER * M_PI * LANDER_SIZE * LANDER_SIZE;
    if(parachute_status==DEPLOYED) { CdA = (DRAG_COEF_LANDER * M_PI + DRAG_COEF_CHUTE * 20) * LANDER_SIZE * LANDER_SIZE; }
    vector3d drag = (- 1.0 / 2)* CdA * atmospheric_density(pos) * v.abs2() * v.norm(); return drag; }

void autopilot ( double Kh, double Kp, double Delta)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
    stabilized_attitude = true;
    // Adjust throttle-->
    double e, h = position.abs()-MARS_RADIUS, v_rad = velocity * position.norm();
    file << simulation_time << ',' << h << ','<<v_rad<<endl;
    e = -(0.5 + Kh*h+ v_rad); //error term
    double P_out = Kp * e;
    if (P_out < -Delta) { throttle = 0; }
    else if(P_out > 1 - Delta) { throttle = 1; }
    else { throttle = Delta + P_out; }

}

void numerical_dynamics(void)
// This is the function that performs the numerical integration to update the
// lander's pose. The time step is delta_t (global variable).
{
    static vector3d previous_position;
    bool verlet = true;
    vector3d drag, grav, thrust, a, next_position;
    double mass;
    mass = UNLOADED_LANDER_MASS + FUEL_DENSITY * FUEL_CAPACITY * fuel;
    grav = -GRAVITY * MARS_MASS * position.norm() / (position.abs2());
    drag = drag_force(velocity, position);
    thrust = thrust_wrt_world();
    a = grav + (thrust + drag) / mass;
    // Integrator select.
    if (verlet==true) {
        if (simulation_time == 0.0) { next_position = position + velocity * delta_t; velocity =velocity + (a * delta_t);}
        else { next_position = (2 * position) - previous_position + (a * delta_t * delta_t);
               velocity = (next_position - position) / delta_t;}
        previous_position = position;
        position = next_position;}
    else{ velocity = velocity + a * delta_t;
          next_position = position + velocity * delta_t;
          position = next_position;}
    // Here we can apply an autopilot to adjust the thrust, parachute and attitude
    //autopilot ( double Kh, double Kp, double Delta)
    if (autopilot_enabled && scenario == 1) autopilot(0.01,69.3,0.6);
    else if (autopilot_enabled && scenario == 3) autopilot(0.00518, 69.3, 0.6);
    // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
    if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "geosychronous orbit";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";

  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = true;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = true;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = true;
    break;

  case 6:
    //geospatial orbit
      position = vector3d( pow(GRAVITY * MARS_MASS * MARS_DAY * MARS_DAY /(4 * M_PI* M_PI), 1.0/ 3.0 ), 0.0, 0.0);
      velocity = vector3d(0.0, sqrt(GRAVITY * MARS_MASS / position.abs()), 0.0);
      orientation = vector3d(0.0, 90.0, 0.0);
      delta_t = 0.1;
      parachute_status = NOT_DEPLOYED;
      stabilized_attitude = false;
      autopilot_enabled = false;
      break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;
  }
}
