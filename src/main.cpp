#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include <math.h>
#include "ukf.h"

using namespace std;

// for convenience
using json = nlohmann::json;
double distance_difference_sum  = 0.0;
double heading_to_target_ud     = 0.0;
double heading_difference_ud    = 0.0;
uint32_t   hunter_hearbeat          = 0;
// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("]");
  if (found_null != std::string::npos) {
    return "";
  }
  else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main()
{
  
  

  uWS::Hub h;

  // Create a UKF instance
  UKF ukf;
  
  double target_x = 0.0;
  double target_y = 0.0;

  h.onMessage([&ukf,&target_x,&target_y](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event

    if (length && length > 2 && data[0] == '4' && data[1] == '2')
    {

      auto s = hasData(std::string(data));
      if (s != "") {
      	
      	
        auto j = json::parse(s);
        std::string event = j[0].get<std::string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          double hunter_x = std::stod(j[1]["hunter_x"].get<std::string>());
          double hunter_y = std::stod(j[1]["hunter_y"].get<std::string>());
          double hunter_heading = std::stod(j[1]["hunter_heading"].get<std::string>());

          string lidar_measurment = j [1]["lidar_measurement"];

          MeasurementPackage meas_package_L;
          istringstream iss_L(lidar_measurment);
    	  long long timestamp_L;

    	  // reads first element from the current line
    	  string sensor_type_L;
    	  iss_L >> sensor_type_L;

          
      	  // read measurements at this timestamp
      	  meas_package_L.sensor_type_ = MeasurementPackage::LASER;
          meas_package_L.raw_measurements_ = VectorXd(2);
          float px;
      	  float py;
          iss_L >> px;
          iss_L >> py;
          meas_package_L.raw_measurements_ << px, py;
          iss_L >> timestamp_L;
          meas_package_L.timestamp_ = timestamp_L;
          

    	  ukf.ProcessMeasurement(meas_package_L);
		 
    	  string radar_measurment = j[1]["radar_measurement"];
          
          MeasurementPackage meas_package_R;
          istringstream iss_R(radar_measurment);
    	  long long timestamp_R;

    	  // reads first element from the current line
    	  string sensor_type_R;
    	  iss_R >> sensor_type_R;

      	  // read measurements at this timestamp
      	  meas_package_R.sensor_type_ = MeasurementPackage::RADAR;
          meas_package_R.raw_measurements_ = VectorXd(3);
          float ro;
      	  float theta;
      	  float ro_dot;
          iss_R >> ro;
          iss_R >> theta;
          iss_R >> ro_dot;
          meas_package_R.raw_measurements_ << ro,theta, ro_dot;
          iss_R >> timestamp_R;
          meas_package_R.timestamp_ = timestamp_R;
          
    	  ukf.ProcessMeasurement(meas_package_R);

	  target_x = ukf.x_[0];
	  target_y = ukf.x_[1];

        // extract state values for read-ability
        double px_p, py_p;
        double p_x  = ukf.x_[0];
        double p_y  = ukf.x_[1];
        double v    = ukf.x_[2];
        double yaw  = ukf.x_[3];
        double yawd = ukf.x_[4];
        // distance from target
        double distance_difference = sqrt((target_y - hunter_y)*(target_y - hunter_y) + (target_x - hunter_x)*(target_x - hunter_x));
        // extrapolate delta time with gain = 0.175
        distance_difference_sum = distance_difference_sum*0.5 + distance_difference;
        double delta_t = distance_difference*0.175 + distance_difference_sum*0.000005;

        // noise free process model
        if (fabs(yawd) > 0.001) 
        {
            px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else 
        {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }
        // set target
        target_x = px_p;
        target_y = py_p;
/*
        hunter_x = target_x;
        hunter_y = target_y;
        hunter_heading = atan2(target_y, target_x);
         */ 
    	  double heading_to_target = atan2(target_y - hunter_y, target_x - hunter_x);
        if(fabs(heading_to_target) < 1000)
        {
    	     while (heading_to_target > M_PI) heading_to_target-=2.*M_PI; 
    	     while (heading_to_target <-M_PI) heading_to_target+=2.*M_PI;
           heading_to_target_ud = heading_to_target;
        }
        else
        {
           heading_to_target = heading_to_target_ud;
        }
    	  //turn towards the target
    	  double heading_difference = heading_to_target - hunter_heading;
        if(fabs(heading_difference) < 1000)
        {
    	     while (heading_difference > M_PI) heading_difference-=2.*M_PI; 
    	     while (heading_difference <-M_PI) heading_difference+=2.*M_PI;
           heading_difference_ud = heading_difference;
        }
        else
        {
           heading_difference = heading_difference_ud;
        }
        cout<<"Hunter heartbeat: "<< hunter_hearbeat++ <<endl;

          json msgJson;
          msgJson["turn"] = heading_difference;
          msgJson["dist"] = distance_difference; 
          auto msg = "42[\"move_hunter\"," + msgJson.dump() + "]";
          // std::cout << msg << std::endl;
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
	  
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }

  });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1)
    {
      res->end(s.data(), s.length());
    }
    else
    {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port))
  {
    std::cout << "Listening to port " << port << std::endl;
  }
  else
  {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}























































































