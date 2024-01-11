/*
 * 
 * Project 4 - Solar Tracker
 * Name - Jayden Wong
 * Student ID - 300618572
 * 
 * */
#include <cmath>
#include "image_pr4.h" 
#include <thread>
#include <chrono>


int redCount = 0; //variable to count red pixels
int greenCount = 0; // counts green pixels
int blueCount = 0; //counts blue pixels
int redPixelTotal = 0; //total count of the red pixels
bool dayCheck = false; //checks if it's daytime or nighttime

int sunStartX, sunStartY, sunEndX, sunEndY; // the x and y values of top and bottom of the sun
//matrices for edge detection
int Gx[3][3] = {{-1, 0, 1},   {-2, 0, 2}, {-1, 0, 1}}; // Gradient matrix for horizontal edge detection
int Gy[3][3] = {{-1, -2, -1}, {0, 0, 0},  {1, 2, 1}}; // Gradient matrix for vertical edge detection

using namespace std;

// Struct to manage and store data related to the path of an orbit
struct Orbit {
	// logged position and time
	std::vector<int> x;  
	std::vector<int> y;
	std::vector<int> t;
	int xc,yc,r;  // center and radius
	int x_sunrise,y_sunrise;
	double omega = 0.1;
} orbit;



// Function to compute the determinant of a 3x3 matrix
double findDeterminant(double a, double b, double c, double d, double e, double f, double g, double h, double i){
	// Calculate the determinant using the formula for a 3x3 matrix
	double determinant = a*(e*i - f*h) - b*(d*i - f*g) + c*(d*h - e*g);
	// Return the calculated determinant value
	return determinant;
}



// determines circle parameters based on data points, used to start predicting the sun
void calculateOrbit() {
    int n = orbit.x.size(); // Number of data points

    // Calculate coefficients for the matrix
    double a1_1 = n;
    double a1_2 = 0;
    double a1_3 = 0;
    for (int i = 0; i < n; i++) {
        a1_3 += cos(orbit.omega * orbit.t[i]); // Sum of cosines of time values
    }

    double a2_1 = 0;
    double a2_2 = n;
    double a2_3 = 0;
    for (int i = 0; i < n; i++) {
        a2_3 += sin(orbit.omega * orbit.t[i]); // Sum of sines of time values
    }

    double a3_1 = 0;
    for (int i = 0; i < n; i++) {
        a3_1 += cos(orbit.omega * orbit.t[i]); // Sum of cosines of time values
    }
    double a3_2 = 0;
    for (int i = 0; i < n; i++) {
        a3_2 += sin(orbit.omega * orbit.t[i]); // Sum of sines of time values
    }
    double a3_3 = n;

    // Calculate constants for the matrix
    double b1 = 0;
    double b2 = 0;
    for (int i = 0; i < n; i++) {
        b1 += orbit.x[i]; // Sum of x-coordinates of data points
        b2 += orbit.y[i]; // Sum of y-coordinates of data points
    }

    double b3_x = 0;
    double b3_y = 0;
    for (int i = 0; i < n; i++) {
        b3_x += orbit.x[i] * cos(orbit.omega * orbit.t[i]); // Sum of products of x-coordinate and cosines of time values
        b3_y += orbit.y[i] * sin(orbit.omega * orbit.t[i]); // Sum of products of y-coordinate and sines of time values
    }
    double b3 = b3_x + b3_y; // Combined sum of the two products

    // Calculate determinants of the matrices
    double determinant = findDeterminant(a1_1, a1_2, a1_3, a2_1, a2_2, a2_3, a3_1, a3_2, a3_3); // Determinant of the original matrix
    double determinantXc = findDeterminant(b1, a1_2, a1_3, b2, a2_2, a2_3, b3, a3_2, a3_3); // Determinant with the first column changed
    double determinantYc = findDeterminant(a1_1, b1, a1_3, a2_1, b2, a2_3, a3_1, b3, a3_3); // Determinant with the second column changed
    double determinantR = findDeterminant(a1_1, a1_2, b1, a2_1, a2_2, b2, a3_1, a3_2, b3); // Determinant with the third column changed

    // Calculate the center coordinates and radius of the circle
    orbit.xc = (determinantXc / determinant); // x-coordinate of the center
    orbit.yc = (determinantYc / determinant); // y-coordinate of the center
    orbit.r = (determinantR / determinant); // radius
}



//CORE & COMPLETION CODE
void locateSun(int time){
	//stores the height and width of the image into different variables
    int height = image.height;
    int width = image.width;
    // for this part, I had help from "red ruby" project from ENGR101 to help me detect the red pixels
	for (int row = 0; row < height; row++){
		for (int col = 0; col < width; col++){
			redCount = int(get_pixel(image, row, col, 0)); //count of red in pixel
		    greenCount = int(get_pixel(image, row, col, 1)); //count of green in pixel
		    blueCount = int(get_pixel(image, row, col, 2)); //count of blue in pixel
		  //if the "redness" is more than the other colours, add it to the total red pixel count.
			if ((redCount > 2 * greenCount) && (redCount > 2 * blueCount)){ 
				redPixelTotal = redPixelTotal + 1; 
			    if(redPixelTotal == 1){ 
					sunStartX = col;
				    sunStartY = row;
				}
				else {
					sunEndX = col;
					sunEndY = row;
				}
			}
		}
	} 
	// if the red pixels are above a certain threshold, that means the sun is on the screen. otherwise it is nighttime.
   if(redPixelTotal > 100){ dayCheck = true; } // the sun is on the screen
   else { dayCheck = false; } //else it's night time
 
   if(dayCheck){
 	   int xCenter = sunStartX + (sunEndX - sunStartX)/2;
	   int yCenter = sunStartY + (sunEndY - sunStartY)/2;
	   orbit.x.push_back(xCenter);
	   orbit.y.push_back(yCenter);  
	   orbit.t.push_back(time);
   }
   //resets redpixeltotal to 0.
   redPixelTotal = 0;
}


//CHALLENGE CODE
// function to detect edges and determine the position of the sun based on the detected edges
void locateEdges(int time) {
    int height = image.height; // Height of the image
    int width = image.width; // Width of the image

    vector<int> edgePointsX; // Stores all the x values of the edges
    vector<int> edgePointsY; // Stores all the y values of the edges

    // Iterate through each pixel in the image to detect edges
    for (int row = 1; row < height - 1; row++) {
        for (int col = 1; col < width - 1; col++) {
            // Determine if the pixel is red
            int redCount = int(get_pixel(image, row, col, 0)); // Red component of the pixel
            int greenCount = int(get_pixel(image, row, col, 1)); // Green component of the pixel
            int blueCount = int(get_pixel(image, row, col, 2)); // Blue component of the pixel

            // Check if the pixel is a red pixel and find the edges
            if ((redCount > 1.4 * greenCount) && (redCount > 1.4 * blueCount)) {
                int GxSum = 0; // Sum of the Gx values
                int GySum = 0; // Sum of the Gy values
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        GxSum = GxSum + get_pixel(image, (row - 1) + i, (col - 1) + j, 0) * Gx[i][j]; // Calculate the sum of Gx values
                        GySum = GySum + get_pixel(image, (row - 1) + i, (col - 1) + j, 0) * Gy[i][j]; // Calculate the sum of Gy values
                    }
                }

                // Check if the pixel is an edge based on the threshold and store the edge point
                double totalSum = sqrt(GxSum * GxSum + GySum * GySum); // Calculate the total sum
                if (totalSum > 80) {
                    edgePointsX.push_back(col); // Store the x coordinate of an edge point
                    edgePointsY.push_back(row); // Store the y coordinate of an edge point
                }
            }
        }
    }

    // Create a counting board for all the pixels
    vector<vector<int>> acc(height, vector<int>(width, 0)); // The accumulator array
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            acc[i][j] = 0; // Initialize all elements to 0
        }
    }

    // Find the position of the sun based on the detected edges
    for (unsigned int i = 0; i < edgePointsX.size(); i++) {
        for (unsigned int j = i; j < edgePointsX.size(); j++) {
            for (unsigned int k = j; k < edgePointsX.size(); k++) {
                // Check if the differences in x and y coordinates are within a certain threshold
                if (abs(edgePointsX[i] - edgePointsX[j]) > 30 && abs(edgePointsX[i] - edgePointsX[k]) > 30 && abs(edgePointsX[j] - edgePointsX[k]) > 30) {
                    continue; // Skip if the threshold is exceeded
                }

                if (abs(edgePointsY[i] - edgePointsY[j]) > 30 && abs(edgePointsY[i] - edgePointsY[k]) > 30 && abs(edgePointsY[j] - edgePointsY[k]) > 30) {
                    continue; // Skip if the threshold is exceeded
                }

                // Compute the position of the potential center of the sun
                double a = 2 * edgePointsX[i] - 2 * edgePointsX[j];
                double b = 2 * edgePointsY[i] - 2 * edgePointsY[j];
                double c = pow(edgePointsX[j], 2) - pow(edgePointsX[i], 2) + pow(edgePointsY[j], 2) - pow(edgePointsY[i], 2);

                double d = 2 * edgePointsX[i] - 2 * edgePointsX[k];
                double e = 2 * edgePointsY[i] - 2 * edgePointsY[k];
                double f = pow(edgePointsX[k], 2) - pow(edgePointsX[i], 2) + pow(edgePointsY[k], 2) - pow(edgePointsY[i], 2);

                if (a * e - b * d != 0) {
                    double xc = (-1 * c * e + f * b) / (a * e - b * d); // Potential center x-coordinate
                    double yc = (-1 * a * f + c * d) / (a * e - b * d); // Potential center y-coordinate
                    if (int(yc) < 900 && int(yc) >= 0 && int(xc) < 900 && int(xc) >= 0) {
                        acc[int(yc)][int(xc)]++; // Increase the vote count in the accumulator array
                    }
                }
            }
        }
    }

    // Find the pixel that has the most votes to determine the position of the sun
    for (int row = 0; row < height; row++) {
        for (int col = 0; col < width; col++) {
            if (acc[row][col] > 30000) {
                // Log the position into the respective vectors
                orbit.x.push_back(col);
                orbit.y.push_back(row);
                orbit.t.push_back(time);
                break;
            }
        }
    }
}




int main(){        
	std::cout<<"Start..."<<std::endl;
	init(2); //change init(x) to what you need. 0 = core, 1 = completion, 2 = challenge.
    for ( int time = 0 ; time < 950; time++){ // how many "ticks" it will iterate through. change the "950" if you want more or less ticks.
      draw_all(time); 
      std::cout<<"Tick = "<<time<<std::endl;
  
      
      if(orbit.x.size() < 30){
		  //locateSun(time);   // core and completion (comment out if doing challenge)
		  locateEdges(time);	// challenge  (comment out if doing core or completion)
		  
      }
      else{
		  // if there is enough data, it can start predicting the sun
		  calculateOrbit();
		  // x and y variables used to calculate orbit parameters
		  double x = orbit.xc + orbit.r * cos(orbit.omega * time);
		  double y = orbit.yc + orbit.r * sin(orbit.omega * time);
		  // Calculate the squared distance from the y-coordinate of the orbit center to the bottom edge of the image. Used to find x-coordinate of sunrise
		  double power = pow((image.height - orbit.yc), 2);
		  // x-coordinate of the sunrise position
		  double xSunrise = orbit.xc - sqrt(abs(orbit.r - power)); 
		  // orbit sunrise positions
		  orbit.x_sunrise = xSunrise; 
		  orbit.y_sunrise = image.height;

		  
		  if(y > image.height){ dayCheck = false; } // if the y value of the orbit parameters is greater than the image height, it is not day. Otherwise, it is day.
		  else { dayCheck = true; } 
		  if(dayCheck){ // Adjust the position of the solar panel based on the position of the sun
			  // the sun is on the screen and so we calculate x and y
			  int xPanel = 0; // x coord of the panel
			  int yPanel = 0; // y coord of the panel
			  get_aim(xPanel, yPanel); //get the current coords of the panel
			  double angle = atan2(y - yPanel, x - xPanel); // Calculate the angle between the current position of the solar panel and the sun
			  move_aim(angle); //moves to the new coords.

		  }
		  else if(!dayCheck){ // Adjust the position of the solar panel based on the sunrise position
			  int xPanel = 0;
			  int yPanel = 0;
			  get_aim(xPanel, yPanel);
			  double angle = atan2(orbit.y_sunrise - yPanel, orbit.x_sunrise - xPanel); // Calculate the angle between the current position of the solar panel and the predicted sunrise position
			  move_aim(angle);
		  }
	  }
      std::this_thread::sleep_for(std::chrono::milliseconds(500)); // Makes program thread to pause for a value of time. You can change this value.
   }

    return 0;
}

