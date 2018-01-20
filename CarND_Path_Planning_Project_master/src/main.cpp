#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"


using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

bool GetBestLane(vector<vector<double>> sensor_fusion, double car_x, double car_y,  double car_s,  double car_d, double car_yaw, double car_speed, int& car_lane, int prev_size,
    vector<double> previous_path_x, vector<double> previous_path_y, vector<double> map_waypoints_s, vector<double> map_waypoints_x, vector<double> map_waypoints_y, double ref_vel) {
  
  vector<bool> bestLane = {true, true, true};

  // std::cout << "NUMBER_VEHICLES: " << sensor_fusion.size() << std::endl;

  // go through the sensor fusion data
  if(sensor_fusion.size() > 0) {

    // std::cout << "In sensor fusion size: " << sensor_fusion.size() << std::endl;

    for(int i=0; i < sensor_fusion.size(); i++) {
      float d = sensor_fusion[i][6];
      int traffic_lane = (int)(d / 4);

      // std::cout << "In lane change Other car d: " << d << "Lane: " << traffic_lane << std::endl;
      // std::cout << "In lane change our car lane: " << car_lane << std::endl;

      // if the car is not in my lane
      if(traffic_lane != car_lane) {
        double vx = sensor_fusion[i][3];
        double vy = sensor_fusion[i][4];
        double traffic_speed = sqrt(vx*vx + vy *vy);
        double other_car_s = sensor_fusion[i][5];

        
        // project s value out in the future
        other_car_s += ((double) prev_size * 0.02 * traffic_speed);
        double predicted_car_s =  car_s + ((double) prev_size * 0.02 * traffic_speed);

        std::cout << "In lane logic predicted s: " << predicted_car_s << std::endl;
        std::cout << "In lane logic other_car_s: " << other_car_s << std::endl;

        if(abs(other_car_s - predicted_car_s) < 45) {

          // std::cout << "Lane not good:" << traffic_lane << std::endl;
          // std::cout << "Diff:" << other_car_s - predicted_car_s << std::endl;
          bestLane[traffic_lane] = false;
        } // end checking other_car_s - predicted_car_s
      } // end checking other car not in my lane
    } // end for loop sensor fusion 
  } // end  sensor fusion size check

  // assign optimum lane

  // for(int i=0; i < bestLane.size(); i++) {
    // std::cout << "Best lane vector:" << i << ": Value" << bestLane[i] << std::endl;
  // }
  
  // std::cout << "Best lane size:" << bestLane.size() << std::endl;

  for(int i=0; i < bestLane.size(); i++) 
  {
    if(bestLane[i] && i != car_lane && abs(i-car_lane) < 2) 
    {
      car_lane = i;
      // std::cout << "Car Lane new:" << car_lane << std::endl;
      break;
    }
  }
} //end function

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }


  
  //1. Define lane and ref velocity
  //start in lane 1
  // Lane 0 is the far left lane, Lane 1 is the middle lane and so on…
  int lane = 1;

  // have a reference velocity 
  double ref_vct = 0.0; // mph


  
  // pass ref velocity and lane variables here
  h.onMessage([&ref_vct, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            json msgJson;

            vector<double> next_x_vals;
            vector<double> next_y_vals;


            // 2. prev path size
            int prev_size = previous_path_x.size();

            // std::cout << "prev size " << prev_size << std::endl;

            // 4. code for sensor fusion data

            // check if the car is in our lane, how close, if too close what to do

            // change the car's s to previous path's last point s
            if (prev_size > 0)
            {
              
                car_s = end_path_s;
            }

            bool too_close = false;

            // cout << "Car Lane: " << lane << endl;

            // Use sensor fusion which is a vector of doubles
            // 2d vector with Car ID and then car x, y position, x, y velocity, car s,d coordinates
            
            // check all cars
            if(sensor_fusion.size() > 0) {

              for (int i=0; i< sensor_fusion.size(); i++ )
              {
              
                  float d = sensor_fusion[i][6];

                  // check if car in same lane as our car
                  if ( (d < (2 + 4*lane + 2)) && (d > (2 + 4*lane - 2)))

                  {

                    // if car in same lane, check speed of car
                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double traffic_speed = sqrt(vx*vx + vy*vy);

                    // get the s value of that car
                    double traffic_s = sensor_fusion[i][5];

                    // now project the s value of the car outwards in time, where will the car be in the future
                    traffic_s+=((double)prev_size*.02*traffic_speed);

                    // if the car  is in front of us and the gap (where car in the future and where we will be in the future) is less than 30 m
                    if ((traffic_s > car_s) && ((traffic_s - car_s) < 30))
                    {
                    
                      too_close = true;

                      // cout << "Car too close, distance:" << traffic_s - car_s << endl;

                      // cout << "Sensor Fusion Size:" << sensor_fusion.size() << endl;

                      // cout << "Start Get best lane****************:"  << endl;

                      GetBestLane(sensor_fusion, car_x, car_y, car_s, car_d, car_yaw, car_speed, lane, prev_size, previous_path_x, previous_path_y, map_waypoints_s, map_waypoints_x, map_waypoints_y, ref_vct);

                      // cout << "too Close =:"  << too_close << endl;
                      // cout << "Car lane final decision =:"  << lane << endl;
                      // cout << "END Get best lane ****************:"  << endl;

                    } // end loop car too close

                  } // end if for checking the d value

              } // end for loop for sensor fusion data

            } // end loop checking if sensor_fusion size > 0


            if (too_close)
            {
              // if the car is too close, decrease the reference velocity by a constant
              ref_vct -= 0.3;
            }
            else if (ref_vct < 49.5)
            {
              // if not too close, gradually increase the ref velocity
              ref_vct+= 0.224;
            }


          

            // END *****code for checking sensor fusion data

            // 3. Code for interpolating waypoints and pushing values to next x vals/y vals
            // create a list of widely spaced (x, y) waypoints evenly spaced at 30 m
            // later interpolate waypoints with a spline

            vector<double> ptsx;
            vector<double> ptsy;

            // reference x, y , yaw states
            // either the start point where the car is or the prev paths end point
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            if (prev_size < 2)
            {
              // use two points that make the path tangent to the car

              // std::cout << "prev size less than 2" << std::endl;
              double prev_car_x = car_x - cos(car_yaw);
              double prev_car_y = car_y - sin(car_yaw);

              ptsx.push_back(prev_car_x);
              ptsx.push_back(car_x);

              ptsy.push_back(prev_car_y);
              ptsy.push_back(car_y);


            }
            // use prev path's end point as starting point
            else
            {
              
            // std::cout << "prev size greater than 2" << std::endl;
            // reference state is prevous path end point
            ref_x = previous_path_x[prev_size - 1];
            ref_y = previous_path_y[prev_size - 1];

            double ref_x_p = previous_path_x[prev_size - 2];
            double ref_y_p = previous_path_y[prev_size - 2];

            //angle the car was heading in using kast two points
            ref_yaw = atan2(ref_y-ref_y_p, ref_x-ref_x_p);


            // use the last and second to last point and angle the car was travelling in
            ptsx.push_back(ref_x_p);
            ptsx.push_back(ref_x);

            ptsy.push_back(ref_y_p);
            ptsy.push_back(ref_y);


            } // end if prev size < 2

            // add 30 m evenly spaced three points 
            // lane is included so that of there is a lane change, the spline is mapped correctly
            vector <double> next_p0 = getXY(car_s +30, (2+4*lane), map_waypoints_s,map_waypoints_x ,map_waypoints_y );
            vector <double> next_p1 = getXY(car_s +60, (2+4*lane), map_waypoints_s,map_waypoints_x ,map_waypoints_y );
            vector <double> next_p2 = getXY(car_s +90, (2+4*lane), map_waypoints_s,map_waypoints_x ,map_waypoints_y );

            // push the points to ptsx and ptsy
            ptsx.push_back (next_p0[0]);
            ptsx.push_back (next_p1[0]);
            ptsx.push_back (next_p2[0]);

            ptsy.push_back (next_p0[1]);
            ptsy.push_back (next_p1[1]);
            ptsy.push_back (next_p2[1]);

            // std::cout << "end created spline points" << std::endl;

            // Make sure the last point of the previous path is at 0, 0 and angle is at 0 degrees
            // transform to local car coordinates

            for (int i=0; i < ptsx.size(); i++)
            {
              
                // shift car reference angle
                double shift_x = ptsx[i] - ref_x;
                double shift_y = ptsy[i] - ref_y;

                ptsx[i] = (shift_x*cos(0-ref_yaw)-shift_y*sin(0-ref_yaw));                
                ptsy[i] = (shift_x*sin(0-ref_yaw)+shift_y*cos(0-ref_yaw));
                

            }

            // std::cout << "end rotate car" << std::endl;

            // at this point the vectors of points have 5 points

            // create a spline
            tk::spline s;

            // set x, y points to the spline
            s.set_points(ptsx, ptsy);


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            // Start defining path


            // First send the remaining previous path points from the last time
            for(int i=0; i < previous_path_x.size(); i++)
            {
              
              next_x_vals.push_back(previous_path_x[i]);
              next_y_vals.push_back(previous_path_y[i]);

            }

            // there are two sets of points 1. The waypoints that make up the spline and 2.path planning points 

            // Calculate x, y points at 30 m and calculate the distance
            double tgt_x = 30.0;
            double tgt_y = s(tgt_x); // get the y value using the spline
            double tgt_dist = sqrt((tgt_x * tgt_x) +  (tgt_y * tgt_y));
                           
            double x_add_on = 0;

            // fill rest of the path planner, always keep 50 points, depends on how many previous path points are there
            for(int i=1; i <= 50-previous_path_x.size(); i++)
            {
              
                // get the x, y points along the path
                double N = tgt_dist/(0.02*ref_vct/2.24);
                double x_point = x_add_on + (tgt_x)/N;
                double y_point = s(x_point);

                x_add_on = x_point;

                double x_ref = x_point;
                double y_ref = y_point;

                // rotate back to normal
                x_point = (x_ref * cos(ref_yaw)-y_ref*sin(ref_yaw));                
                y_point = (x_ref *sin(ref_yaw) +y_ref*cos(ref_yaw));
             
                // std::cout << "rotate back to normal completed" << std::endl;

                x_point += ref_x;
                y_point += ref_y;

                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }

            // std::cout << "end adding all points" << std::endl;


            // END defining path

            
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
