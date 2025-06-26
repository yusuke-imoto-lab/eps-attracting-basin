Column             | Description 
-------------------|------------------------------------------------------------------------
InitializationTime | Date-Time of the start of the forecast (YYYY-mm-dd HH:MM:SS UTC)
ForecastTime       | Date-Time of the forecast (YYYY-mm-dd HH:MM:SS UTC)
Horizon            | Forecast horizon in hour since the start (from 0 to +39)
Member             | Id of the ensemble member (0 to 20)
type               | Type of ensemble member (unperturbed or perturbed)      
lon                | Longitude of the center in decimal degrees
lat                | Latitude of the center in decimal degrees
X                  | Longitude of the center in meters (projection: JGD2011 / UTM zone 53N)
Y                  | Latitude of the center in meters (projection: JGD2011 / UTM zone 53N)
prmsl              | Pressure reduce to the mean sea level (in Pa)
Xs                 | Scaled longitude (Standard score; unitless)   
Ys                 | Scaled latitude (Standard score; unitless)
X1                 | X - First X position, scaled by the standard deviation at each IT
Y1                 | Y - First Y position, scaled by the standard deviation at each IT
id1                | id1 (YYYYmmddHHMMSS_YYYYmmddHHMMSS_em)
id2                | id2 (ITid_FTid_em)