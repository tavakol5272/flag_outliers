# Movement Probability Outlier Flagger

MoveApps

Github repository: *github.com/movestore/????*

## Description
This app detects and flags movement outliers in animal tracking data by modelling step length and turning angle probabilities for each track, using a user-defined probability type and threshold.

## Documentation

The app detects movement outliers in animal tracking data based on movement probabilities derived from step length and turning angle.

The user can control how outliers are defined via the `Threshold` and `Probability type` settings:
**Threshold**: Controls how strict outlier detection is (default: 0.05), flagging only the lowest-probability locations as outliers.
**Probability type**: lets the user choose which movement probability to use for outliers (step_turn, delta_step, delta_turn, joint, or custom), based on different combinations of step length and turning angle.

First, the app splits the data by track ID and detects outliers separately for each track. 
Then it applies several pre-cleaning steps: it removes empty geometries, missing timestamps, 
 time lags shorter than 60 seconds, extremely high speeds (top 1%), and all locations from the first and last day. 
 After cleaning, if an individual track has fewer than three locations, outlier detection for that track is skipped.
 
The app computes step length between consecutive locations and turning angle between successive steps, 
then estimates a 2D distribution over step length × turning angle and 1D distributions for changes in step length and turning angle. 
It builds a 2D histogram (raster) of step length vs. turning angle and usesthe selected probability measures (joint, step_turn, delta_step, delta_turn, custom) to flag outliers. 

For each track, two plots are produced: All locations (coloured by probability) and Kept vs removed (showing which points are flagged as outliers).

Before returning, the function prints a small table of outlier locations.

On return, the app outputs the original data as a move2 object with additional columns:
step_length_mv, turning_angle_mv, step_turn_prob, delta_step_prob, delta_turn_prob, joint_prob, custom_prob, log_prob, outlier_percentile, is_outlier, and is_na_prob. 


### Application scope
#### Generality of App usability
This App was developed for any taxonomic group. 

#### Required data properties
The App was developed and tested with Sigfox tracking data, but it works with any kind of move2 object.

### Input type
`move2::move2_loc`

### Output type
`move2::move2_loc`

### Artefacts
`Plots` :
For each track, two plots are produced:
1.	**All locations (probability)**: Points colored by log probability.
2.	**Kept vs removed** : outliers in red and Kept points in blue.

**logs output**: Before returning, the function prints a small table of outlier locations.

**Returned data** : returns the input `move2` object with additional columns for movement metrics, probabilities, and outlier flags.

### Settings

`Threshold:`: Probability percentile as a fraction (0–1). e.g.: 0.05 = flag lowest 5% as outliers.

`Probability type`: Dropdown to select select which movement probability is used for outlier detection. Options:
**step_turn**: Probability of the current step length + turning angle combination.(Default)
**joint**: Combined probability of step_turn × delta_step × delta_turn (most restrictive)
**delta_step**: Probability of the change in step length from one step to the next.
**delta_turn**: Probability of the change in turning angle from one step to the next.
**custom**: Combined probability of step_turn × delta_step only.

### Changes in output data

The App returns the input `move2` object with additional columns:
`step_length_mv`, `turning_angle_mv`, `step_turn_prob`, `delta_step_prob`, `delta_turn_prob`, `joint_prob`, `custom_prob`, `log_prob`, `outlier_percentile`, `is_outlier`, `is_na_prob`. 

### Most common errors
Track id(s) with `NA`: Outlier detection requires valid track IDs for all locations.  
Too few locations after cleaning: If a track has fewer than 3 locations after pre-cleaning, outlier detection for that track is skipped.


### Null or error handling
If the input is `NULL` or has 0 rows, the App returns `NULL` and logs an info message.  
Tracks with fewer than 3 locations after cleaning are returned unchanged (no outlier flags)
