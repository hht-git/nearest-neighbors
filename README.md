The purpose of the R code is to find K nearest neighbors for the spatial points with longitude and latitude coordinates.

<p align="center">
  <b></b><br>
  <img src="Cities and Nearest Stations.gif">
</p>

* Two orders of magnitude faster than st_nn function in *nngeo* package.
* Longitude/latitude coordinates to xyz coordinates transformation.
* Chord to arc and arc to chord convertion on earth great circle.
* 3 dimensional k-d tree searching.
* Distance returned same as the result by Haversine function in *geosphere* package.
* Methods to use non-vectorized function in *dplyr*'s mutate.


