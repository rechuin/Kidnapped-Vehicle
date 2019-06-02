/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */


/*
#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <sstream>

#include "helper_functions.h"
*/


#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

//using std::string;
//using std::vector;
//using std::normal_distribution;
using namespace std;

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
  default_random_engine gen;  // initialize random engine
  normal_distribution<double> dist_x(x,std[0]);  
  normal_distribution<double> dist_y(y,std[1]);
  normal_distribution<double> dist_theta(theta,std[2]);
  int i; 
  Particle current_particle;
  for(i=0; i<num_particles; i++)
  {
    // random x, y, theta refer to their normal distrubution
    current_particle.x = dist_x(gen);  // refer to lesson 6-6
    current_particle.y = dist_y(gen);
    current_particle.theta = dist_theta(gen);
    current_particle.weight = 1.0;
    
    // associate each elements into corrensponds vector
    particles.push_back(current_particle);
    weights.push_back(current_particle.weight);
  }
  
  // flag setting
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
  default_random_engine gen;  // initialize random engine
  int i;
  
  for (i = 0; i < num_particles; i++)
  {
    double pred_x;
    double pred_y;
    double pred_theta;
    
    // prediction step: need to calculate each x,y based on yaw angle theta whether equal to 0
    // refer to lessons 6-8&9
    if(fabs(yaw_rate) < 0.0001)
    {
      pred_x = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
      pred_y = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
      pred_theta = particles[i].theta;
    }
    else
    {
      pred_x = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      pred_y = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      pred_theta = particles[i].theta + yaw_rate*delta_t;
    }
    
    // random variation for each elements
    normal_distribution<double> dist_x(pred_x,std_pos[0]);  
    normal_distribution<double> dist_y(pred_y,std_pos[1]);
    normal_distribution<double> dist_theta(pred_theta,std_pos[2]);
    
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations, double sensor_range) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  int i,j;
  
  // each observation will have a nearst landmark ID.
  
  //std::cout << "observation number is : " << observations.size() << std::endl;
  //std::cout << "prediction land mark number is : " << predicted.size() << std::endl;
  
  
  for (i = 0; i<observations.size(); i++)
  {
    double nearst_dist = sensor_range*sqrt(2); // Maximum distance can be squre root of 2 times
    int nearst_landmark_id = -1;  // initialize a nearst landmark id
    
    for (j=0; j<predicted.size(); j++)  // cycle prediction point for matching nearst distance with observation
    {
      double current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);  // calculate the distance between observation and prediction
      
      // compare with the neast distance, if so give the corresponds id to nearst landmark id
      if(current_dist < nearst_dist) 
      {
        nearst_dist = current_dist;
        nearst_landmark_id = predicted[j].id;
      }
    }
    // update the observation id using landmark id
    observations[i].id = nearst_landmark_id;
    //std::cout << "observation id is : " << observations[i].id << std::endl;
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
  
  int i,j;
  default_random_engine gen;  // initialize random engine
  double weight = 0.0;
  
  for(i=0; i<num_particles; i++)
  {
    // restore partical attributes
    double particle_x = particles[i].x;
    double particle_y = particles[i].y;
    double particle_theta = particles[i].theta;
    particles[i].weight = 1.0;
    /* */
    /* step1: transformate observation coordinate to map coordinate */
    /* */
    
    // creat tranformatd observation vector
    vector<LandmarkObs> transformed_observations;
    
    for(j=0; j<observations.size(); j++)
    {
      LandmarkObs transf_obs;  // decleration vector
      
      transf_obs.id = j;  // mark observation tranformed ID
      // transfrom observation coordinate to map based on eacb patical 
      transf_obs.x = particle_x + (cos(particle_theta)*observations[j].x) - (sin(particle_theta)*observations[j].y);  
      transf_obs.y = particle_y + (sin(particle_theta)*observations[j].x) - (cos(particle_theta)*observations[j].y);
      transformed_observations.push_back(transf_obs);  // associate into vector
    }
    
    /* */
    /*Step 2: Filter map landmarks to keep only those which are in the sensor_range of current particle. Push them to predictions vector.*/
    /* */
    vector<LandmarkObs> predicted_landmarks;
    
    // ensure all distance between partical and landmark less than to sensor range
    for (j = 0; j < map_landmarks.landmark_list.size(); j++) {
      if ((fabs((particle_x - map_landmarks.landmark_list[j].x_f)) <= sensor_range) && (fabs((particle_y - map_landmarks.landmark_list[j].y_f)) <= sensor_range)) 
      { 
        LandmarkObs pred_lm;
        pred_lm.id = map_landmarks.landmark_list[j].id_i;
        pred_lm.x = map_landmarks.landmark_list[j].x_f;
        pred_lm.y = map_landmarks.landmark_list[j].y_f;
        predicted_landmarks.push_back(pred_lm);
      }
    }
    
    /* */
    /*Step 3: Associate observations to lpredicted andmarks using nearest neighbor algorithm.*/
    /* */
    
    //Associate observations with predicted landmarks
    dataAssociation(predicted_landmarks, transformed_observations, sensor_range);
    
    /* */
    /*Step 4: Calculate the weight of each particle using Multivariate Gaussian distribution.*/
    /* */
    
    double sigma_x = std_landmark[0];
    double sigma_y = std_landmark[1];
    double normalizer;
    double term1;
    double term2;
    double term_exp;
    double total_weight = 0.0;
    double weight_int;
    int k,l;
    
    for (k=0; k<transformed_observations.size(); k++)
    {
      for (l=0; l<predicted_landmarks.size(); l++)
      {
        if (transformed_observations[k].id == predicted_landmarks[l].id)
        {
          normalizer = 1.0/(2.0 * M_PI* sigma_x * sigma_y);
          term1 = pow((transformed_observations[k].x - predicted_landmarks[l].x), 2)/(2.0 * sigma_x * sigma_x);
          term2 = pow((transformed_observations[k].y - predicted_landmarks[l].y), 2)/(2.0 * sigma_y * sigma_y);
          term_exp = exp(-1.0 * (term1 + term2));
          weight_int = normalizer * term_exp; // refer to lessons 6-20
          if(weight_int == 0)
          {
            weight_int = 0.001;
          }
          particles[i].weight *= weight_int;
          /*
          std::cout << "sigma_x is : " << sigma_x << std::endl;
          std::cout << "sigma_y is : " << sigma_y << std::endl;
          std::cout << "transformed_observations x is : " << transformed_observations[k].x << std::endl;
          std::cout << "transformed_observations y is : " << transformed_observations[k].y << std::endl;
          std::cout << "predicted_landmarks x is : " << predicted_landmarks[l].x << std::endl;
          std::cout << "predicted_landmarks y is : " << predicted_landmarks[l].y << std::endl;
          std::cout << "normalizer is : " << normalizer << std::endl;
          std::cout << "term1 is : " << term1 << std::endl;
          std::cout << "term2 is : " << term2 << std::endl;
          std::cout << "term_exp is : " << term_exp << std::endl;
          std::cout << "weight_int is : " << weight_int << std::endl;
          std::cout << "****************************************************"<< std::endl;
          */
        }
      }
    }
    total_weight += particles[i].weight;
    
    /* */
    /*Step 5: Normalize the weights of all particles since resmapling using probabilistic approach.*/
    /* */
    
    for (int i = 0; i < particles.size(); i++) {
    particles[i].weight /= total_weight;
    weights[i] = particles[i].weight;
    }
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> resampled_particles;

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;
	
	//Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
	
	for (int i = 0; i < particles.size(); i++) {
		uniform_real_distribution<double> random_weight(0.0, max_weight_2);
		beta += random_weight(gen);

	  while (beta > weights[current_index]) {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % num_particles;
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
particles = resampled_particles;

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
