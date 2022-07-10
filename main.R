source('functions.R')

library(tictoc)
library(nngeo)
library(geosphere)
library(plotly)
library(htmlwidgets)

# get map of USA and states
usa_sf <- map_data('usa') %>% 
  st_as_sf(coords = c('long', 'lat'), 
           crs = st_crs(4326))

state <- map_data("state")
state_sf <- state %>% 
  st_as_sf(coords = c('long', 'lat'), 
           crs = st_crs(4326))


sf_use_s2(FALSE) # to avoid edge crosses edge error.
usa_bound <- usa_sf %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") %>% 
  st_make_valid() %>%
  st_union()


SEARCH_RADIUS_KM <- 500 # searching radius in km
K <- 10 # return first K nearest neighbors

# create dummy data for cities and weather stations roughly in US

# start a small set of data to be shown on map
N_CITIES <- 20
N_STATIONS <- 500

# create cities in usa
cities_usa_sf <- points_in_usa(N_CITIES, usa_bound, TRUE) %>% 
  mutate(index = row_number(),
         city_name = paste0('city_', index))

# create weather stations in usa
weather_stations_usa_sf <- points_in_usa(N_STATIONS, usa_bound, TRUE) %>% 
  mutate(index = row_number(),
         station_name = paste0('station_', index))

# get K nearest stations and the distances from the cities 
connection <- nn(cities_usa_sf, weather_stations_usa_sf,
                 k = K, 
                 search_radius_km = SEARCH_RADIUS_KM) %>%
  # combine original city/station info
  inner_join(cities_usa_sf, by = c('x.idx' = 'index')) %>% 
  inner_join(weather_stations_usa_sf, by = c('y.idx' = 'index'))

# create a simple map to show the nearest stations connected to the cities
p <- ggplot() + 
  geom_polygon(data=state, 
               aes(x = long, y = lat, group = group, text = region), 
               fill = 'white', color = 'black', size = 0.1) + 
  geom_sf(data = weather_stations_usa_sf,
          aes(text = paste0('candidate station: ', station_name,
                            '\n(', longitude, ', ', latitude, ')')),
          color = 'blue', size = 0.5, alpha = 0.2) +
  geom_segment(data = connection,
               aes(x = longitude.x, y = latitude.x, xend = longitude.y, yend = latitude.y),
               size = 0.1) +
  geom_point(data=connection,
             aes(x=longitude.y, y=latitude.y, 
                 text = paste0('selected station: ', station_name,
                               '\n(', longitude.y, ', ', latitude.y, ')')),
             color = 'blue', size = 0.5) +
  geom_point(data=connection,
             aes(x=longitude.x, y=latitude.x, 
                 text = paste0('city: ', city_name,
                 '\n(', longitude.x, ', ', latitude.x, ')')),
             color = 'red', size = 1) +
  xlab(NULL) + ylab(NULL) +
  ggtitle(paste('Cities and', K, 'Nearest Weather Stations')) +
  theme(plot.title = element_text(hjust = 0.5,size=12))
p
ggsave('Cities and Nearest Stations.jpg',
       units = 'in', width = 9, height = 6, dpi = 300, type = 'cairo')

# interactive map, mouse tooltip shows state name and city/station index
# zoom by mouse scrolling
p_plotly <- ggplotly(p, tooltip = 'text') %>% config(scrollZoom = TRUE)
p_plotly

htmlwidgets::saveWidget(p_plotly, "tmp.html", selfcontained = TRUE)
file.rename("tmp.html", "Cities and Nearest Stations.html")

# try second method for distances
connection2 <- nn_2(cities_usa_sf, weather_stations_usa_sf,
                    k = K, 
                    search_radius_km = SEARCH_RADIUS_KM) %>%
  # combine original city/station info
  inner_join(cities_usa_sf, by = c('x.idx' = 'index')) %>% 
  inner_join(weather_stations_usa_sf, by = c('y.idx' = 'index'))

R_m <- 6371000 # average earth radius in meters

# compare the distances returned by different functions
# the distance functions in geosphere package can not be used directly in dplyr's mutate
# try several ways to call these functions 
# these methods also work for calling other non-vectorized functions
distance_compare <- connection %>% select(everything(), dist_km) %>% 
  inner_join(connection2 %>% 
               rename(dist_km_2 = dist_km)) %>% 
  # by row method
  mutate(distHaversine = by(., 1:nrow(.), function(row) {
    distHaversine(unlist(row$geometry.x), 
                  unlist(row$geometry.y), 
                  R_m)/1000})) %>% 
  # map method
  mutate(distGeo = pmap(select(., starts_with('geometry')), 
                        ~ distGeo(c(unlist(..1)[1], unlist(..1)[2]), 
                                  c(unlist(..2)[1], unlist(..2)[2]), 
                                  R_m/1000))) %>% 
  # rowwise method
  rowwise() %>% 
  mutate(distMeeus  = distMeeus (unlist(geometry.x), 
                                 unlist(geometry.y), 
                                 R_m)/1000) %>% 
  ungroup() %>% 
  # the function takes matrix as arguments
  mutate(distCosine = distCosine(matrix(unlist(geometry.x), ncol=2, byrow =T), 
                                 matrix(unlist(geometry.y), ncol=2, byrow =T), 
                                 R_m)/1000) %>% 
  select(city_name, 
         station_name, 
         dist_km, 
         dist_km_2, 
         distHaversine, 
         distGeo, 
         distMeeus,
         distCosine)


# when datasets are small, e.g. 20 cities and 500 stations, 
# the funcion nn finishes almost instantly

# create large datasets to test the performance of the function
N_CITIES2 <- 1000
N_STATIONS2 <- 100000

# create cities/stations.
cities_usa_sf_2 <- points_in_usa(N_CITIES2, usa_bound) %>% 
  mutate(index = row_number(),
         city_name = paste0('city_', index))

# create weather stations in usa
weather_stations_usa_sf_2 <- points_in_usa(N_STATIONS2, usa_bound) %>% 
  mutate(index = row_number(),
         station_name = paste0('station_', index))

# get K nearest stations and the distances from the cities 
# large dataset, still very fast.
tic('search by function nn')
connection2 <- nn(cities_usa_sf_2, weather_stations_usa_sf_2,
                  k = K, 
                  search_radius_km = SEARCH_RADIUS_KM) %>%
  # combine original city/station info
  inner_join(cities_usa_sf_2, by = c('x.idx' = 'index')) %>% 
  inner_join(weather_stations_usa_sf_2, by = c('y.idx' = 'index'))
toc()
# 2.51 sec elapsed


# st_nn function in nngeo is very slow
tic('search by st_nn in nngeo')
distances_nngeo <- st_nn(cities_usa_sf_2, 
                         weather_stations_usa_sf_2, 
                         k = K, 
                         maxdist = SEARCH_RADIUS_KM * 1000, # maxdist needs unit meter
                         returnDist = T)
toc()
# 351.36 sec elapsed
