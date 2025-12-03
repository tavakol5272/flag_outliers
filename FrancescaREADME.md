# Outlier Detection based on tag error

MoveApps

Github repository: *github.com/movestore/????*

## Description
This app detects and flags movement outliers in animal tracking data based on an expected tag error rate. 
Given a user-defined probability type and threshold, it uses step length, turning angle and related movement probabilities to flag the least likely locations as outliers in each individual.



## Documentation

**Notes:**   
  1-Before running the app, the user must change Attribute defining track ID in the option tab(Configure Movebank Location Datasource) from "Combination of animal and deployment (default)"
to "Animal", because subsequent analyses will be performed by this track ID.

  2-This app is designed for low-frequency tracking data (e.g. Sigfox) with fixes at least 5 minutes apart in the tag and is not suitable for high resolution data.


The app assumes that a certain fraction of locations from a tag may be erroneous (e.g. 5%) and uses movement probabilities to identify those least plausible points per individual.

The user controls how outliers are defined via the `Threshold` and `Probability type` settings:

- **Threshold**: Fraction of locations expected to be erroneous for this tag (default: 0.05). e.g., if the threshold is defined to 0.05 so that about the 5% least probable locations are flagged as outliers. 

- **Probability type**: Chooses which movement probability is used for outlier detection (`Step turn`, `Delta step`, `Delta turn`, `Joint`, or `Custom`), based on different combinations of step length 
  and turning angle.  
  
-**Outlier handling**: "Choose whether outliers are only flagged or removed"(`Flag outliers (keep them in the data)` , `Remove outliers (delete them from the data)`)
  
-**Notes:**:  
  1-If the threshold is not known from the tag engineers, an appropriate probability type and threshold should be chosen by checking visually and trying different settings.  
  2- **Strongly recommend** that, if you want to remove outliers, you first run the app with flagging only. Once you are confident that the settings are appropriate, you can rerun it with outlier removal enabled.

First, the app splits the data by individuals; user should change the Attribute defining track ID in option tab (Configure Movebank Location Datasource) when getting the data from Movebank to "Animal"
(not the default one:combinationof animal and deployment). The subsequent analyses will be performed by this track ID.

Then it applies several pre-cleaning steps. Removes:  
  -empty geometries and timestamps  
  -time lags shorter than 60 seconds: fixes are at least 5 minutes apart in the tag, so points closer than 60 seconds to each other are wrong fixes   
  -extremely high speeds (top 1%):mostly caused by tag/GPS errors rather than true animal movement  
  -all locations from the first and last day: tagging is stressful and the last day data are typically incomplete  
 
After cleaning, if an individual has fewer than three locations, outlier detection for that individual is skipped because turning angle cannot be computed with less than 3 points.
 
The app computes step length between consecutive locations and turning angle between successive steps, 
then estimates a 2D distribution over step length × turning angle and 1D distributions for changes in step length and turning angle. 
It builds a 2D histogram (raster) of step length vs. turning angle and uses the selected probability measures (joint, step_turn, delta_step, delta_turn, custom) to flag outliers. 

For each individual, two plots are produced: All locations (coloured by probability of being or not an outlier) and Kept vs removed (showing which points are flagged as outliers).

Before returning, the function prints a small table of outlier locations.

On return, the app outputs the original data as a move2 object with additional columns:  
(step_length_mv, turning_angle_mv, step_turn_prob, delta_step_prob, delta_turn_prob, joint_prob, custom_prob, log_prob, outlier_percentile, is_outlier, and is_na_prob). 


### Application scope
#### Generality of App usability
This App can be used for any taxonomic group, but it is designed for low-frequency GPS tracking data (e.g. Sigfox GPS tags with fixes every ~5–10 minutes), 
not for high-resolution (e.g. 1 Hz) data or triangulation based.


#### Required data properties
The App was developed and tested with Sigfox tracking data, but it works with any kind of move2 object.

### Input type
`move2::move2_loc`

### Output type
`move2::move2_loc`

### Artefacts
`Plots` :  
For each individual, two plots are produced:  
1.	**All locations (probability)**: Points colored by log probability.  
2.	**Kept vs removed** : outliers in red and Kept points in blue.

**logs output**: Before returning, the function prints a small table of outlier locations.

**Returned data** : returns the input `move2` object with additional columns for movement metrics, probabilities, and outlier flags(showing in column: `is_outlier`= TRUE).

### Settings

`Threshold:`: Probability percentile as a fraction (0–1). e.g.: 0.05 = flag the 5% of locations with the lowest probability to be real locations.

`Probability type`: Dropdown to select which movement probability is used for outlier detection. Options:  
**Step turn**: Probability of the current step length + turning angle combination.(Default)  
**Delta step**: Probability of the change in step length from one step to the next.  
**Delta turn**: Probability of the change in turning angle from one step to the next.  
**Joint**: Combined probability of step_turn × delta_step × delta_turn (most restrictive)  
**Custom**: Combined probability of step_turn × delta_step only.

`Outlier handling`: Radiobuttons to choose whether outliers are only flagged or removed. options:    
-**Flag outliers**: keep them in the data    
-**Remove outliers:** delete them from the data

### Changes in output data

The App returns the input `move2` object with additional columns:

`step_length_mv`: Step length between consecutive locations    
`turning_angle_mv`: Turning angle between successive steps    
`step_turn_prob`: Probability of the observed combination of step length and turning angle   
`delta_step_prob`: Probability of the change in step length from one step to the next   
`delta_turn_prob`: Probability of the change in turning angle from one step to the next   
`joint_prob`: Combined probability of `step_turn_prob × delta_step_prob × delta_turn_prob`   
`custom_prob`: Combined probability of `step_turn_prob × delta_step_prob`  
`log_prob`: Log10-transformed probability used for coloring the plots    
`outlier_percentile`: percentile where higher values correspond to more unlikely locations (i.e. lower movement probability)   
`is_outlier`: the location is an outlier = TRUE  
`is_na_prob`: the selected probability is `NA`= TRUE

### Most common errors
-Track id(s) with `NA`: Outlier detection requires valid track IDs for all locations.  
-Too few locations after cleaning: If an individual has fewer than 3 locations after pre-cleaning, outlier detection for that individual is skipped because turning angle cannot be computed with less than 3 points.


### Null or error handling
If the input is `NULL` or has 0 rows, the App returns `NULL` and logs an info message.  
Individuals with fewer than 3 locations after cleaning are returned unchanged (no outlier flags)
