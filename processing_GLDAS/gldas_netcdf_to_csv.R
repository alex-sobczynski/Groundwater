library(ncdf4) ##unpacking netcdfs
library(rgdal) ##reading shapefile
library(rgeos) ##spdf management
library(sf) ##other stuff
library(tictoc)
#read shapefile (whole country)
afghan <- readOGR("C:/Users/asobc/PycharmProjects/uc-research-2022/data/smoist/scripts/country/AFG_adm0.shp")
#create new CRS and assign to shapefile
CRS.new <-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
proj4string(afghan) <- CRS.new
# set path and filename
ncpath <- "C:/Users/asobc/PycharmProjects/ClimateProject_2023/cygwin_fold_2023/all_netcdfs/"
##ncpath <- "D:/ncdf_fldas_2016/"
##make the dates so the year can be changed
#year <- "2013"
#day <- c("01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28")
#day2 <-c("29","30")
#day3 <-c("01","02","03","04","05","06","07","08","09","10","11","12")
#date <- c(paste(year,"02",day, sep=""),paste(year,"03",day3, sep=""))
##date <- c( "20130216", paste(year,"02",day, sep=""), paste(year,"03",day, sep=""),paste(year,"03",day2, sep=""),paste(year,"03","31", sep=""),paste(year,"04",day, sep=""),paste(year,"04",day2, sep=""),paste(year,"05",day, sep=""),paste(year,"05",day2, sep=""),paste(year,"05","31", sep=""),paste(year,"06","01", sep=""))
#make all the file names
files <- paste0(ncpath,list.files(ncpath))
#test <- "C:/Users/asobc/PycharmProjects/ClimateProject_2023/cygwin_fold_2023/all_netcdfs/GLDAS_CLSM025_DA1_D.A20230101.022.nc4"
#check <- substr(test,102,109)
#finding 
dname <- "Lwnet_tavg"  # note: tmp means temperature (not temporary)
aname <- "Qle_tavg"
bname <- "Qh_tavg"
cname <- "Qg_tavg"
ename <- "Evap_tavg"
fname <- "Qs_tavg"
gname <- "AvgSurfT_tavg"
hname <- "Qsb_tavg"
iname <- "Qsm_tavg"
jname <- "SnowT_tavg"
kname <- "SWE_tavg"
lname <- "SnowDepth_tavg"
mname <- "SoilMoist_S_tavg"
nname <- "SoilMoist_RZ_tavg"
oname <- "TWS_tavg"
pname <- "GWS_tavg"
qname <- "SoilMoist_P_tavg"
rname <- "TVeg_tavg"
#now we loop
for(file in files) {
  # Open the nc file
  t <- try(ncin <- nc_open(file))
  if("try-error" %in% class(t)) {
    # Error occured
    print(paste0("THIS FILE BROKE: ", file))
    next
  }
  print(paste0("Working with file: ", file))
  
  #
  # # open a netCDF file
  # ncin <- nc_open(file)
  nc_close(ncin)
}
print("Done.")

#To do our logic on the files:

for(file in files) {
  # Check if the file already converted
  csvpath <- "C:/Users/asobc/Harris-Public-Policy Dropbox/Alex Sobczynski/WorldBank_AFG_Climate/data/GLDAS/csv_2023/"
  time <-substr(file,102,109)
  csvname <- paste("gldas_",time,".csv",sep="")
  csvfile <- paste(csvpath, csvname, sep="")
  if(file.exists(csvfile)){
    # Skip the file
    next
  }
  
  # Manual skipping of files
  ##if(file == "D:/ncdf_fldas_2020/FLDAS_NOAH001_G_CA_D.A20200415.001.nc"){
  ##  print("WARNING!!! Manually skipping file")
  ##  next
  ##}
  
  # Start the stopwatch
  tic()
  
  # Open the nc file
  t <- try(ncin <- nc_open(file))
  if("try-error" %in% class(t)) {
    # Error occured
    print(paste0("THIS FILE BROKE: ", file))
    next
  }
  print(paste0("Working with file: ", file))
  
  # Get the longitudes and latitude
  lon <- ncvar_get(ncin,"lon")
  lat <- ncvar_get(ncin,"lat")
  
  # Create the arrays
  d_array <- ncvar_get(ncin,dname)
  dlname <- ncatt_get(ncin,dname,"long_name")
  dunits <- ncatt_get(ncin,dname,"units")
  fillvalued <- ncatt_get(ncin,dname,"_FillValue")
  # get a
  a_array <- ncvar_get(ncin,aname)
  alname <- ncatt_get(ncin,aname,"long_name")
  aunits <- ncatt_get(ncin,aname,"units")
  fillvaluea <- ncatt_get(ncin,aname,"_FillValue")
  # get c
  c_array <- ncvar_get(ncin,cname)
  clname <- ncatt_get(ncin,cname,"long_name")
  cunits <- ncatt_get(ncin,cname,"units")
  fillvaluec <- ncatt_get(ncin,cname,"_FillValue")
  #get b
  b_array <- ncvar_get(ncin,bname)
  blname <- ncatt_get(ncin,bname,"long_name")
  bunits <- ncatt_get(ncin,bname,"units")
  fillvalueb <- ncatt_get(ncin,bname,"_FillValue")
  # get e
  e_array <- ncvar_get(ncin,ename)
  elname <- ncatt_get(ncin,ename,"long_name")
  eunits <- ncatt_get(ncin,ename,"units")
  fillvaluee <- ncatt_get(ncin,ename,"_FillValue")
  # get f
  f_array <- ncvar_get(ncin,fname)
  flname <- ncatt_get(ncin,fname,"long_name")
  funits <- ncatt_get(ncin,fname,"units")
  fillvaluef <- ncatt_get(ncin,fname,"_FillValue")
  # get g
  g_array <- ncvar_get(ncin,gname)
  glname <- ncatt_get(ncin,gname,"long_name")
  gunits <- ncatt_get(ncin,gname,"units")
  fillvalueg <- ncatt_get(ncin,gname,"_FillValue")
  # get h
  h_array <- ncvar_get(ncin,hname)
  hlname <- ncatt_get(ncin,hname,"long_name")
  hunits <- ncatt_get(ncin,hname,"units")
  fillvalueh <- ncatt_get(ncin,hname,"_FillValue")
  # get i
  i_array <- ncvar_get(ncin,iname)
  ilname <- ncatt_get(ncin,iname,"long_name")
  iunits <- ncatt_get(ncin,iname,"units")
  fillvaluei <- ncatt_get(ncin,iname,"_FillValue")
  # get j
  j_array <- ncvar_get(ncin,jname)
  jlname <- ncatt_get(ncin,jname,"long_name")
  junits <- ncatt_get(ncin,jname,"units")
  fillvaluej <- ncatt_get(ncin,jname,"_FillValue")
  # get k
  k_array <- ncvar_get(ncin,kname)
  klname <- ncatt_get(ncin,kname,"long_name")
  kunits <- ncatt_get(ncin,kname,"units")
  fillvaluek <- ncatt_get(ncin,kname,"_FillValue")
  # get l
  l_array <- ncvar_get(ncin,lname)
  llname <- ncatt_get(ncin,lname,"long_name")
  lunits <- ncatt_get(ncin,lname,"units")
  fillvaluel <- ncatt_get(ncin,lname,"_FillValue")
  # get m
  m_array <- ncvar_get(ncin,mname)
  mlname <- ncatt_get(ncin,mname,"long_name")
  munits <- ncatt_get(ncin,mname,"units")
  fillvaluem <- ncatt_get(ncin,mname,"_FillValue")
  # get n
  n_array <- ncvar_get(ncin,nname)
  nlname <- ncatt_get(ncin,nname,"long_name")
  nunits <- ncatt_get(ncin,nname,"units")
  fillvaluen <- ncatt_get(ncin,nname,"_FillValue")
  # get o
  o_array <- ncvar_get(ncin,oname)
  olname <- ncatt_get(ncin,oname,"long_name")
  ounits <- ncatt_get(ncin,oname,"units")
  fillvalueo <- ncatt_get(ncin,oname,"_FillValue")
  # get p
  p_array <- ncvar_get(ncin,pname)
  plname <- ncatt_get(ncin,pname,"long_name")
  punits <- ncatt_get(ncin,pname,"units")
  fillvaluep <- ncatt_get(ncin,pname,"_FillValue")
  # get q
  q_array <- ncvar_get(ncin,qname)
  qlname <- ncatt_get(ncin,qname,"long_name")
  qunits <- ncatt_get(ncin,qname,"units")
  fillvalueq <- ncatt_get(ncin,qname,"_FillValue")
  # get r
  r_array <- ncvar_get(ncin,rname)
  rlname <- ncatt_get(ncin,rname,"long_name")
  runits <- ncatt_get(ncin,rname,"units")
  fillvaluer <- ncatt_get(ncin,rname,"_FillValue")
  
  ###These should be rectangular
  # replace netCDF fill values with NA's
  a_array[a_array==fillvaluea$value] <- NA
  b_array[b_array==fillvalueb$value] <- NA
  c_array[c_array==fillvaluec$value] <- NA
  d_array[d_array==fillvalued$value] <- NA
  e_array[e_array==fillvaluee$value] <- NA
  f_array[f_array==fillvaluef$value] <- NA
  g_array[g_array==fillvalueg$value] <- NA
  h_array[h_array==fillvalueh$value] <- NA
  i_array[i_array==fillvaluei$value] <- NA
  j_array[j_array==fillvaluej$value] <- NA
  k_array[k_array==fillvaluek$value] <- NA
  l_array[l_array==fillvaluel$value] <- NA
  m_array[m_array==fillvaluem$value] <- NA
  n_array[n_array==fillvaluen$value] <- NA
  o_array[o_array==fillvalueo$value] <- NA
  p_array[p_array==fillvaluep$value] <- NA
  q_array[q_array==fillvalueq$value] <- NA
  r_array[r_array==fillvaluer$value] <- NA
  # create dataframe -- reshape data
  # matrix (nlon*nlat rows by 2 cols) of lons and lats
  lonlat <- as.matrix(expand.grid(lon,lat))
  # vectors of values
  a_vec <- as.vector(a_array)
  b_vec <- as.vector(b_array)
  c_vec <- as.vector(c_array)
  d_vec <- as.vector(d_array)
  e_vec <- as.vector(e_array)
  f_vec <- as.vector(f_array)
  g_vec <- as.vector(g_array)
  h_vec <- as.vector(h_array)
  i_vec <- as.vector(i_array)
  j_vec <- as.vector(j_array)
  k_vec <- as.vector(k_array)
  l_vec <- as.vector(l_array)
  m_vec <- as.vector(m_array)
  n_vec <- as.vector(n_array)
  o_vec <- as.vector(o_array)
  p_vec <- as.vector(p_array)
  q_vec <- as.vector(q_array)
  r_vec <- as.vector(r_array)
  
  remove(a_array)
  remove(b_array)
  remove(c_array)
  remove(d_array)
  remove(e_array)
  remove(f_array)
  remove(g_array)
  remove(h_array)
  remove(i_array)
  remove(j_array)
  remove(k_array)
  remove(l_array)
  remove(m_array)
  remove(n_array)
  remove(o_array)
  remove(p_array)
  remove(q_array)
  remove(r_array)
  
  # create dataframe and add names
  df01 <- data.frame(cbind(lonlat,a_vec,b_vec,c_vec,d_vec,e_vec,f_vec,g_vec,h_vec,i_vec,j_vec,k_vec,l_vec,m_vec,n_vec,o_vec,p_vec,q_vec,r_vec))
  names(df01) <- c("lon","lat",paste0(aname,"_",time),paste0(bname,"_",time),paste0(cname,"_",time),paste0(dname,"_",time),paste0(ename,"_",time),paste0(fname,"_",time),paste0(gname,"_",time),paste0(hname,"_",time),paste0(iname,"_",time),paste0(jname,"_",time),paste0(kname,"_",time),paste0(lname,"_",time),paste0(mname,"_",time),paste0(nname,"_",time),paste0(oname,"_",time),paste0(pname,"_",time),paste0(qname,"_",time),paste0(rname,"_",time))
  #cropping in afghanistan
  df02 <- subset(df01, df01$lon > 60 & df01$lon < 75 & df01$lat > 29 & df01$lat < 39)
  ##create an spdf for df02
  spdf <-SpatialPointsDataFrame(coords =df02[, c("lon","lat")], data=df02,
                                proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  remove(df02)
  spdf2 <- spdf[!is.na(over(spdf, as(afghan, "SpatialPolygons"))), ]
  df03 <- as.data.frame(spdf2)
  df04 <- subset(df03, select = -c(lat.1,lon.1))
  
  # clean up file from memory
  nc_close(ncin)
  remove(file)
  remove(df01)
  remove(df02)
  remove(df03)
  print(names(df04))  
  # Move the csv file
  write.csv(df04,csvfile)
  print(paste0("Created file: ", csvfile))
  
  
  toc()
}
print("Done.")
