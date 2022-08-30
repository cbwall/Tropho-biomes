# temp loggers script

#clear work space
rm(list=ls())
ls()

# load packages
library(reshape2)
library(plyr)
library(lubridate)
library(tools)
library(dplyr)

# remove all unwanted rows and columns
# split date and time column

##### grab files in a list
files <- list.files(path="data/hobo", pattern = "csv$", full.names = T)
files

##### what are the file names, sans extensions
file.names<-file_path_sans_ext(list.files(path="data/hobo/", pattern = "csv$", full.names = F))
file.names


############ formatting all data in for loop
for(i in 1:length(files))
{
  data<-read.csv(files[i], sep=",", header=FALSE)
  df<-data[c(-1:-2),c(-1,-4:-7)] # removes "trash" columns
  #df<-data[, c(2:3)] # reorganize
  colnames(df)<-c("Date.time", "Raw.Temp")
  df$Raw.Temp<-as.numeric(as.character(df$Raw.Temp))
  df$timestamp<-mdy_hms(as.character(df$Date.time)) # corrects date format
  df$timestamp<-strptime(df$timestamp, format="%Y-%m-%d %H:%M:%S")
  df<-df[, c(3,2)] # remove old date-time column
  df<-df[!(df$timestamp < "2022-06-25 00:00:00"),] # start at this time
  df<-df[!(df$timestamp > "2022-08-17 00:00:00"),] # end at this time
  # make.names(assign(paste("Yos2022_",file.names[i], sep=""), df)) if you want a pattern on names
  make.names(assign(paste(file.names[i], sep=""), df))
  # makes each df[i] as dataframe with specific file-name
  # write.csv(df.out, file=paste("trim",file.names[i])) # makes .csvs for output
}
# this is the end of the loop

########## see the files you've made, as a list, then grab those that fill your pattern 
ls()
files<-as.data.frame(mget(ls(pattern = "SN.*"))) # with SN as patterns in files from for loops
names(files) # see number of columns, and what these are

data_index<-c(1,(seq(2,42,2))) # these are the columns we will want: timestamp + raw data **change '42' to number of columns in your dataframe, specifying here to select 'every other column'

temps<-as.data.frame(c(files[, data_index])) # here is the data we want, now in df alone
#strip name to SN only, the "\\..*" = keep everything before the period but not after, the "*"
# why? "." is considered as a special character, so must use a double backslash in front (i.e. \\).
names(temps) = gsub(pattern = "\\..*", replacement = "", x = names(temps)) 

colnames(temps)[1]="timestamp" # rename the single column for time

### apply data callibration to each column
# ex: temp$calib.SN10209515<-(temp$SN10209515 * "value") # Temp logger x: SN x

write.csv(temps, "output/Yos_2022_hobodata.csv")



#### #### #### #### #### #### #### #### #### 
#### calibrate data
#### #### #### #### #### #### #### #### #### 

#clear work space
rm(list=ls())
ls()

####  pull in the periods when the loggers were being calibrated
files <- list.files(path="data/hobo", pattern = "csv$", full.names = T)

##### what are the file names, sans extensions
file.names<-file_path_sans_ext(list.files(path="data/hobo/", pattern = "csv$", full.names = F))

############ formatting all data in for loop
for(i in 1:length(files))
{
  data<-read.csv(files[i], sep=",", header=FALSE)
  df<-data[c(-1:-2),c(-1,-4:-7)] # removes "trash" columns
  #df<-data[, c(2:3)] # reorganize
  colnames(df)<-c("Date.time", "Raw.Temp")
  df$Raw.Temp<-as.numeric(as.character(df$Raw.Temp))
  df$timestamp<-mdy_hms(as.character(df$Date.time)) # corrects date format
  df$timestamp<-strptime(df$timestamp, format="%Y-%m-%d %H:%M:%S")
  df<-df[, c(3,2)] # remove old date-time column
  df<-df[!(df$timestamp < "2022-08-19 10:30:00"),] # start at this time
  df<-df[!(df$timestamp > "2022-08-19 13:00:00"),] # end at this time
  make.names(assign(paste(file.names[i], sep=""), df))
}
# this is the end of the loop

########## see the files you've made, as a list, then grab those that fill your pattern 
ls()
files<-as.data.frame(mget(ls(pattern = "SN.*")))
names(files) 
data_index<-c(1,(seq(2,42,2)))
log.cal<-as.data.frame(c(files[, data_index]))
names(log.cal) = gsub(pattern = "\\..*", replacement = "", x = names(log.cal)) 

colnames(log.cal)[1]="timestamp" # rename the single column for time

###### "log.cal" is the file with all the hobo data

# bring in the YSI data
YSI.cal<-read.csv("data/hobo/YSI calibration/calibration.csv")

colnames(YSI.cal)<-c("Date", "Time", "YSI.Temp")
YSI.cal$Date<-as.Date(YSI.cal$Date, format="%m/%d/%y")
YSI.cal$timestamp<-as.POSIXct(paste(YSI.cal$Date, YSI.cal$Time), format="%Y-%m-%d %H:%M:%S")

# select the columns we want
YSI.cal<-YSI.cal %>%
  select(timestamp, YSI.Temp)


##### merge the dataframes by timestamp
calib.df<-merge(YSI.cal, log.cal ,by="timestamp")

##### ##### ##### ##### ##### ##### ##### 
# specify dataframe and calibration file
df<- calib.df  # set dataframe
y<-calib.df$YSI.Temp # set "y" <<< this is the calibrated value or the "standard"

#####  
# run regressions with equation of lines, export figure

#  *note* pdf() and dev.off() outside the loop will generate the cumulative figure
pdf(paste("output/Yos2022_hobo_calibr_fig.pdf",sep=""))
par(mfrow=c(2,2))

mods<-list()
for (i in 3:23) { # runs through all columns of loggers
  x<-df[,i] 
  m <- lm(y~ x, data = df) # run the regression
  plot(y ~ x, ylab="YSI Temp C", xlab="Logger Temp C", main=colnames(df)[i])
  abline(m)
  eq <-substitute(italic(y) == a + b~italic(x)*","~~italic(r)^2~"="~r2,
                  list(a = format(coef(m)[1], digits = 4),
                       b = format(coef(m)[2], digits = 4),
                       r2 = format(summary(m)$r.squared, digits = 3)))
  legend("topleft", legend=eq, bty="n")
  # to save model output with SN... below gives SN, intercept and slope
  mods[[i]]<-c(colnames(df)[i], format(coef(m)[1], digits = 4), format(coef(m)[2], digits = 4))
}
dev.off()

# models run from above, exported into a dataframe
models.df<-do.call(rbind.data.frame, mods)
colnames(models.df)<-c("SN", "Intercept", "slope")
write.csv(models.df, "output/Yos_2022_hobo_calib_equat.csv")

#### #### #### #### #### BOOOOMMMMMMM! #### #### #### #### #### ####

# apply calibrations

Temp.df<-read.csv("output/Yos_2022_hobodata.csv")
Temp.df$timestamp<-as.POSIXct(Temp.df$timestamp, format="%Y-%m-%d %H:%M:%S") 

# new dataframe for calibrated data
Temp.cal<-as.data.frame(Temp.df$timestamp); colnames(Temp.cal)<-"timestamp"


Temp.cal$SN10339184<-Temp.df$SN10339184*1.036 + -0.4593 


