/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles

  // Set the size of the particles and weights vectors 
  particles.resize(num_particles);
  weights.resize(num_particles);

  std::default_random_engine gen;
  //Create normal(gaussian) distributions for x, y, theta
  std::normal_distribution<double> dist_x(x, std[0]);
  std::normal_distribution<double> dist_y(y, std[1]);
  std::normal_distribution<double> dist_theta(theta, std[2]);

  // Create each individual particle and add it to the particles vector
  // set the all weights in the weights vector to 1
  for(int i = 0; i < this->num_particles; ++ i){

    Particle part;
    part.id = i;
    part.x = dist_x(gen);
    part.y = dist_y(gen);
    part.theta = dist_theta(gen);
    part.weight = 1;

    particles[i] = part;
    weights[i] = 1;
  }

  is_initialized = true; 
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // Generate random Gaussian noise to each particle
  std::default_random_engine gen;
  std::normal_distribution<double> dist_x(0.0, std_pos[0]);
  std::normal_distribution<double> dist_y(0.0, std_pos[1]);
  std::normal_distribution<double> dist_theta(0.0, std_pos[2]);

  // for everty part in particle
  for(std::vector<Particle>::iterator part = particles.begin(); part != particles.end(); part ++){

    double x_pred, y_pred, theta_pred;
    

    if(fabs(yaw_rate) < 0.0000001){
      x_pred = (*part).x + velocity * delta_t * cos((*part).theta);
      y_pred = (*part).y + velocity * delta_t * sin((*part).theta);
      theta_pred = (*part).theta;
    }
    else{
      x_pred = (*part).x + velocity / yaw_rate * (sin((*part).theta + yaw_rate * delta_t)-sin((*part).theta));
      y_pred = (*part).y + velocity / yaw_rate * (cos((*part).theta)-cos((*part).theta + yaw_rate * delta_t));
      theta_pred = (*part).theta + yaw_rate * delta_t;
    }

    // add random Gaussian noise.
    x_pred += dist_x(gen);
    y_pred += dist_y(gen);
    theta_pred += dist_theta(gen);

    // update particle location with prediction
    (*part).x = x_pred;
    (*part).y = y_pred;
    (*part).theta = theta_pred;
  }
}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */

  // For each observation Find the closest prediction 

  // for each obsercation in the observations vector
  for(auto obs = observations.begin(); obs != observations.end(); obs ++){
    double nearest_dist = std::numeric_limits<double>::max();
    
    // for every prediction in the predicted vector
    for(auto pred = predicted.begin(); pred != predicted.end(); pred ++){
      double distance = dist((*obs).x, (*obs).y, (*pred).x, (*pred).y);
      if(distance < nearest_dist){
        nearest_dist = distance;
        (*obs).id = (*pred).id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  int i = 0;
  for(std::vector<Particle>::iterator part = particles.begin(); part != particles.end(); part ++){
    
    std::vector<LandmarkObs> observations_transformed;
    // Transform observation vector
    // For each obsevation change it's perspective
    for(auto obs = observations.begin(); obs != observations.end(); obs ++){
      LandmarkObs temp;
      temp.x = (*part).x + cos((*part).theta) * (*obs).x - sin((*part).theta) * (*obs).y;
      temp.y = (*part).y + sin((*part).theta) * (*obs).x + cos((*part).theta) * (*obs).y;
      observations_transformed.push_back(temp);
    }

    // Make prediction vector
    // find landmarks within the sensor range only
    std::vector<LandmarkObs> predicted;
    for(auto landmark = map_landmarks.landmark_list.begin(); landmark != map_landmarks.landmark_list.end(); landmark ++){

      double distance = dist((*landmark).x_f, (*landmark).y_f, (*part).x, (*part).y);

      if(distance <= sensor_range){
        LandmarkObs pred;
        pred.id = (*landmark).id_i;
        pred.x = (*landmark).x_f;
        pred.y = (*landmark).y_f;
        predicted.push_back(pred);
      }
    }
    // Assosiate data
    dataAssociation(predicted, observations_transformed);

    double new_w = 1.0;

    for(auto obs = observations_transformed.begin(); obs != observations_transformed.end(); obs ++){

      double mu_x = 0;
      double mu_y = 0;


      for(auto pred = predicted.begin(); pred != predicted.end(); pred ++){

        if((*pred).id == (*obs).id){
          mu_x = (*pred).x;
          mu_y = (*pred).y;
        }
      }

      // product
      new_w *= multivariable_gauss(std_landmark[0], std_landmark[1], (*obs).x, (*obs).y, mu_x, mu_y);
    }

    (*part).weight = new_w;
    weights[i++] = new_w;
  }
}



void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // find the largest weight.
  double largest_w = 0.0;

  for(auto part = weights.begin(); part != weights.end(); part++){

    if(*part > largest_w){
      largest_w = *part;
    }    
  }
  std::default_random_engine random(time(NULL));
  std::uniform_int_distribution<int> uniform_dist_int(0, num_particles);
  std::uniform_real_distribution<double> uniform_dist_double(0, 2*largest_w);

  double beta = 0.0;
  int index = uniform_dist_int(random);

  std::vector<Particle> particles_resample;
  std::vector<double> weights_resample;


  for(int i=0; i<num_particles; i++)
  {
    beta += uniform_dist_double(random);
    while(weights[index] < beta)
    {
      beta -= weights[index];
      index = (index+1) % num_particles;
    }
    Particle part_resample;
    part_resample.id = i;
    part_resample.x = particles[index].x;
    part_resample.y = particles[index].y;
    part_resample.theta = particles[index].theta;
    part_resample.weight = particles[index].weight;

    particles_resample.push_back(part_resample);
    weights_resample.push_back(part_resample.weight);
  }

  particles = particles_resample;
  weights = weights_resample;
}


void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}