# Script with functions to calculate the distance between two points on the earth's
# surface using the Haversine formula.

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
haversine <- function(long1, lat1, long2, lat2) {
  
  # Ensure Lats and Longs are in radians
  long1 <- deg2rad(long1)
  lat1 <- deg2rad(lat1)
  long2 <- deg2rad(long2)
  lat2 <- deg2rad(lat2)
  
  # Calculate geodesic distance based on havesine formala
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1 - a))
  d = R * c
  return(d) # Distance in km
}
