files:

okada.py- calc_deformation(alpha, strike, depth, dip,  strike_width, dip_width, dislocation, x, y)
uses okada wrapper to calculate deformation at x,y from a fault patch

plotting.py- creates a latitude and longitude meshgrid and scatter plots

read_model.py- reads in slip files and creates patch objects

run_all.py- uses surface grid and patch objects to loop through every patch and every obs point and plots

to do:
1) plot patches for visualization
2) check okada inputs
3) loop through patches and points in run_all and create superposition deformation map
4) plot differences of gaus and wang
5) plot given stations for both models on and offshore
