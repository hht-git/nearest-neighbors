library(tidyverse)
library(sf)
library(RANN)

# return closest point index for objects x and y 
# size of x should be equal or smaller than size of y
# x, y can be either dataframes with longitude/latitude columns or sf object of points
# if the column names contain 'lon'/'x' or 'lat'/'y' in dataframe (case insensitive),
# latitude/longitude column index will be auto searched, 
# search_radius > 5000km will increase running time dramatically
# 1000 rows in x, 100000 rows in y, k = 10, search_radius = 500km
# st_nn 300+ seconds, this function around 3 seconds

nn <- function(x, y, k = 1, search_radius_km = 1000,
                              xlon_col = 0, xlat_col = 0, 
                              ylon_col = 0, ylat_col = 0,
                              R = 6371) {
  # R = 6378.137 # radius of the Earth in km
  # R = 6371 # average radius of the Earth in km
  
  ctoa <- function(c) {return (R * asin(c / 2 / R) * 2)} # chord to arc on earth's great circle
  atoc <- function(a) {return (R * sin(a / 2 / R) * 2)} # arc to chord on earth's great circle
  
  # transform longitude/latitude to xyz coordinates
  calc_xyz <- function(df, lon_col = 0, lat_col = 0) {
    # print(class(df))
    if (class(df)[1] == 'sf') {
      df <- st_coordinates(df) %>% as.data.frame()
      lon_col = 1
      lat_col = 2
    }
    # print(lon_col)
    # print(lat_col)
    if (lon_col  ==  0) lon_col <- match(TRUE, str_detect(tolower(names(df)), 'lon'))
    if (lat_col  ==  0) lat_col <- match(TRUE , str_detect(tolower(names(df)), 'lat'))
    if (is.na(lon_col)) lon_col <- match('x', tolower(names(df)))
    if (is.na(lat_col)) lat_col <- match('y', tolower(names(df)))
    if (is.na(lon_col)) stop('column longitude or x is missing')
    if (is.na(lat_col)) stop('column latitude or y is missing')
    df %>% mutate(x = R * cos(.[[lon_col]] * pi / 180) * cos(.[[lat_col]] * pi / 180), 
                  y = R * sin(.[[lon_col]] * pi / 180) * cos(.[[lat_col]] * pi / 180), 
                  z = R * sin(.[[lat_col]] * pi / 180)) %>% 
      dplyr::select(x, y, z)
  }
  
  nn2(calc_xyz(y, lon_col = ylon_col, lat_col = ylat_col), 
      calc_xyz(x, lon_col = xlon_col, lat_col = xlat_col), 
      k = k, 
      searchtype = 'radius', radius = atoc(search_radius_km)) %>% 
    as.data.frame() %>% 
    mutate(x.idx = row_number()) %>% 
    pivot_longer(!x.idx, names_to = c('.value'), names_pattern = 'nn.(idx|dist).*[0-9]*') %>% 
    filter(idx !=  0) %>% 
    mutate(dist = ctoa(dist)) %>% 
    dplyr::select(x.idx, y.idx = idx, dist_km = dist)
}


# second way by reprojecting the coordinates to equal distance projection
nn_2 <- function(x, y, k = 1, search_radius_km = 1000,
                                xlon_col = 0, xlat_col = 0, 
                                ylon_col = 0, ylat_col = 0) {
  # transform longitude/latitude to 2d meters by projection to crs 5070
  calc_meter <- function(df, lon_col = 0, lat_col = 0) {
    if (class(df)[1] != 'sf') {
      if (lon_col  ==  0) lon_col <- match(TRUE, str_detect(tolower(names(df)), 'lon'))
      if (lat_col  ==  0) lat_col <- match(TRUE , str_detect(tolower(names(df)), 'lat'))
      if (is.na(lon_col)) lon_col <- match('x', tolower(names(df)))
      if (is.na(lat_col)) lat_col <- match('y', tolower(names(df)))
      if (is.na(lon_col)) stop('column longitude or x is missing')
      if (is.na(lat_col)) stop('column latitude or y is missing')
      sf <- df %>% st_as_sf(coords = c(lon_col,lat_col), crs = st_crs(4326))
    } else sf <- df
    sf %>% st_transform(5070) %>% st_coordinates(df) %>% as.data.frame() %>% rename(x=X,y=Y)
  }
  nn2(calc_meter(y, lon_col = ylon_col, lat_col = ylat_col), 
      calc_meter(x, lon_col = xlon_col, lat_col = xlat_col), 
      k = k, treetype = "kd",
      searchtype = 'radius', search_radius_km * 1000) %>% 
    as.data.frame() %>% 
    mutate(x.idx = row_number()) %>% 
    pivot_longer(!x.idx, names_to = c('.value'), names_pattern = 'nn.(idx|dist).*[0-9]*') %>% 
    filter(idx !=  0) %>% 
    mutate(dist = dist / 1000) %>% 
    dplyr::select(x.idx, y.idx = idx, dist_km = dist)
}

# create sf points in usa
# suppress the message "st_intersection assumes that they are planar"
points_in_usa <- function(n_points, usa_bound, exactly_in = FALSE) {
  # exactly_in = TRUE, the returned points are in usa boundary
  # exactly_in = FALSE, the returned points are in usa boundary box (rectangle)
  usa_bound_box <- usa_bound %>% st_bbox()
  if (exactly_in) {
  tibble(index = 1:(n_points * 2), 
         longitude = round(runif(n_points * 2, 
                                 min = usa_bound_box$xmin, 
                                 max = usa_bound_box$xmax), 
                           3), 
         latitude = round(runif(n_points * 2, 
                                min = usa_bound_box$ymin, 
                                max = usa_bound_box$ymax), 
                          3)) %>% 
    st_as_sf(coords = c('longitude', 'latitude'), 
             remove = FALSE,
             crs = st_crs(4326)) %>% 
    st_intersection(usa_bound) %>% 
    sample_n(n_points) %>%
    suppressWarnings() %>%
    suppressMessages()
  } else {
    tibble(index = 1:n_points, 
           longitude = round(runif(n_points, 
                                   min = usa_bound_box$xmin, 
                                   max = usa_bound_box$xmax), 
                             3), 
           latitude = round(runif(n_points, 
                                  min = usa_bound_box$ymin, 
                                  max = usa_bound_box$ymax), 
                            3)) %>% 
      st_as_sf(coords = c('longitude', 'latitude'), 
               remove = FALSE,
               crs = st_crs(4326))
  }
}

