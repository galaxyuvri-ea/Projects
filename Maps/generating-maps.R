# R script to generate alpha diversity plots for the tick-viruses dataset
# Authors: Alfred Ssekagiri, Phionah Tushabe
# Link: https://github.com/galaxyuvri-ea/Projects/blob/master/Maps/generating-maps.R
# Created: 05th-October-2022
# Modified: 05th-October-2022

# clear environment
rm(list = ls())

# load required packages
library(sf)
library(ggplot2)

# download shape file from https://data.unhcr.org/en/documents/details/83043
# download sample data from https://simplemaps.com/data/ug-cities

# THIS IS THE ONLY PART WE NEED TO CHANGE

# specify PATH to required file; first the shape file
shape_file="~/Desktop/MISC/Uganda_Districts-2020---136-wgs84/Uganda_Districts-2020---136-wgs84.shp"
# secondly specify PATH to csv data for what you want to display, 
# it should have atleast three columns; lat, lng, city -- change as appropriate
sample_data_file="~/Downloads/ug.csv"

# NO NEED TO CHANGE ANYTHING AFTER THIS POINT

# import the shape file that has been downloaded from link above
uga<-st_read(shape_file)
# import the sample data
df<-read.csv(sample_data_file)
df<-df[1:10,] # comment this out, we used it to test the script

# set up a map object using the information in shae file imported above
map<-ggplot(data = uga) + geom_sf(fill="white",color="#E5E4E2") + theme_void()
# add sample data to the map
map<-map+geom_point(data = df,aes(x=lng,y=lat),color="blue")+
  geom_text(data=df, aes(x=lng,y=lat, label=city))+
  theme(legend.position = "none")
map