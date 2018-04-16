/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	num_particles = 50;
	
	default_random_engine generator;
	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> distribution_x(x, std[0]);
	
	// TODO: Create normal distributions for y and theta.
	normal_distribution<double> distribution_y(y, std[1]);
	normal_distribution<double> distribution_theta(theta, std[2]);
	double weight = 1.0;
	Particle p;
	for(int i=0; i<num_particles; i++){
		double sample_x, sample_y, sample_theta;
		
		// TODO: Sample  and from these normal distrubtions like this: 
		// sample_x = dist_x(gen);
		// where "gen" is the random engine initialized earlier.
		sample_x = distribution_x(generator);
		sample_y = distribution_y(generator);
		sample_theta = distribution_theta(generator);
		
		weights.push_back(weight);
		
		p.id = i;
		p.x = sample_x;
		p.y = sample_y;
		p.theta = sample_theta;
		p.weight = weights.at(i);
		
		particles.push_back(p);
	}
	
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	// x_f = x_0 + v/yaw_rate ( sin(theta_0 + yaw_rate*dt) - sin(theta_0))
	// y_f = y_0 + v/yaw_rate ( cos(theta_0) - cos(theta_0 + yaw_rate*dt))
	// theta_f = theta_0 + yaw_rate*dt
	
	default_random_engine generator;
	// This line creates a normal (Gaussian) distributions.
	normal_distribution<double> distribution_x(0, std_pos[0]);
	normal_distribution<double> distribution_y(0, std_pos[1]);
	normal_distribution<double> distribution_theta(0, std_pos[2]);
	Particle p;
	for(int i=0; i<num_particles; i++){
	
		p = particles.at(i);
		
		if(fabs(yaw_rate) <= 0.0001){
			p.x += velocity * delta_t * cos(p.theta);
			p.y += velocity * delta_t * sin(p.theta);
		} else {
			p.x += (( velocity / yaw_rate) * ( sin( p.theta + ( delta_t*yaw_rate)) - sin(p.theta)));
			p.y += (( velocity / yaw_rate) * ( cos(p.theta) - cos( p.theta + ( delta_t*yaw_rate))));
			p.theta += delta_t*yaw_rate;
		}
		
		p.x += distribution_x(generator);
		p.y += distribution_y(generator);
		p.theta += distribution_theta(generator);
		
		particles[i] = p;
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for(int i=0; i< observations.size();i++){
		LandmarkObs obs = observations[i];
		
		double minimum_dist = numeric_limits<double>::max();
		
		int landmark_id = -1;
		
		for(int j=0;j<predicted.size(); j++){
			LandmarkObs prd = predicted[j];
			
			double current_dist = dist(obs.x, obs.y, prd.x, prd.y);
			
			if (current_dist < minimum_dist){
				minimum_dist = current_dist;
				landmark_id = prd.id;
			}
		}
		
		observations[i].id = landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	for(int i = 0; i < num_particles; i++){

		vector<LandmarkObs> predictions;

		for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
			if(fabs(map_landmarks.landmark_list[j].x_f - particles[i].x) <= sensor_range && fabs(map_landmarks.landmark_list[j].y_f - particles[i].y) <= sensor_range){
				predictions.push_back(LandmarkObs{ map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f});
			}
		}

		vector<LandmarkObs> transformed_obs;
		for(int j = 0; j < observations.size(); j++){
			double trans_x = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
			double trans_y = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
			transformed_obs.push_back(LandmarkObs{ observations[j].id, trans_x, trans_y });
		}
		
		dataAssociation(predictions, transformed_obs);
		
		particles[i].weight = 1.0;

		for(int j = 0; j < transformed_obs.size(); j++){
			double obs_x, obs_y, prd_x, prd_y;
			obs_x = transformed_obs[j].x;
			obs_y = transformed_obs[j].y;

			int associated_prediction = transformed_obs[j].id;

			for(int k = 0; k < predictions.size(); k++){
				if(predictions[k].id == associated_prediction){
					prd_x = predictions[k].x;
					prd_y = predictions[k].y;
				}
			}

			double obs_w = ( 1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp( -( pow(prd_x-obs_x,2)/(2*pow(std_landmark[0], 2)) + (pow(prd_y-obs_y,2)/(2*pow(std_landmark[1], 2)))));

			particles[i].weight *= obs_w;
			weights[i] = particles[i].weight;
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	default_random_engine generator;
	uniform_real_distribution<double> uniform(0, 1);
	
	vector<Particle> resampled;
	double max_weight = 0.0;
	for(int i=0;i<num_particles;i++){
		if(weights[i] > max_weight){
			max_weight = weights[i];
		}
	}
	
	int index = (int)(uniform(generator) * num_particles);
	double beta = 0.0;
	for(int i=0;i<num_particles;i++){
		beta += uniform(generator) * 2 * max_weight;
		
		while(beta > weights[index]){
			beta -= weights[index];
			index = (index + 1) % num_particles;
		}
		
		resampled.push_back(particles[index]);
	}
	
	particles = resampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
