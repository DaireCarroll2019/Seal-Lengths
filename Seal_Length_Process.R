#Daire Carroll, Gothenburg University, 2021
#a pipeline based on the curved_length function to process a folder full of shape file subfolders
#users will be required to make a final descision for seals > 1.6 m

########################################

my.dir = paste("") #here paste your file location eg: ~/Desktop/Drones/All_Shape_files

setwd(my.dir)
files = dir()

Overwrite = TRUE

processed_data = matrix(ncol = 8)
colnames(processed_data) = c("Polygon","Length_A","Width_A","lat","lon","Length_M","Width_M","File")

for(i in 1:length(files)){ 
  
  folder = files[i]
  setwd(paste(my.dir,"/",folder,sep = ""))
  
  print(folder)
  
  files2 = dir()
  
  for(k in 1:length(files2 )){
    if(endsWith(files2[k],".shp")==TRUE){
      shape =  files2[k]
    }
  }
  
  print(shape)
  
  seals = st_read(shape)
  seals = st_transform(seals,23032)
  
  seals_dimensions = matrix(ncol = 5, nrow = length(st_geometry(seals)))
  seals_dimensions[,1] = c(1:length(st_geometry(seals)))
  
  for(j in 1:length(st_geometry(seals))){
    pol1 = seals[j,]
    lens = curved_length(pol1, plt = FALSE)
    seals_dimensions[j,2:ncol(seals_dimensions)] = lens
  }

  joined = cbind(seals,seals_dimensions[,1],seals_dimensions[,2])
  
  names(joined)[names(joined) == "seals_dimensions...1."] = "Length_A"
  names(joined)[names(joined) == "seals_dimensions...2."] = "Width_A"

  if(Overwrite == TRUE){
    st_write(joined, shape, delete_layer = TRUE)
  }
  
  if(length(which(colnames(pol1) == "Length_M")) == 1 && length(which(colnames(pol1) == "Width") == 1 )){
    seals_dimensions = cbind(seals_dimensions,seals$Length_M,seals$Width,rep(folder,length(seals_dimensions[,1])))
  }else{
    seals_dimensions = cbind(seals_dimensions,rep(NA, length(seals_dimensions[,1])),rep(NA, length(seals_dimensions[,1])),rep(folder,length(seals_dimensions[,1])))
  }
  
  setwd(my.dir)

  processed_data = rbind(processed_data, seals_dimensions)
   
}

abnormal = which(processed_data[,2] > 1.6)

if(length(abnormal)>=1){
  for(i in 1:length(abnormal)){
    print(processed_data[abnormal[i],2])
    folder = processed_data[abnormal[i],ncol(processed_data)]
    setwd(paste("~/Desktop/Drones/All_Shape_files","/",folder,sep = ""))
 
    files2 = dir()
    
    for(k in 1:length(files2 )){
      if(endsWith(files2[k],".shp")==TRUE){
        shape =  files2[k]
      }
    }
    
    seals = st_read(shape)
    seals = st_transform(seals,23032)

    pol = seals[processed_data[abnormal[i],1],]
    curved_length(pol, plt = TRUE)
    
    Quest = readline("Does this polygon measument look reasonable... Y to leave in dataset, N to remove")
    if(regexpr(Quest, 'y', ignore.case = TRUE) == 1){
      continue = TRUE
    }else if (regexpr(Quest, 'n', ignore.case = TRUE) == 1){
      print("Removing automatic measurments from data set, consider splitting the polygon?")
      
      processed_data[abnormal[i],2:3] = NA
      
      seals[processed_data[abnormal,1],7:8] = NA
    
      Quest = readline("Overwrite shape file? Y/N")
      if(regexpr(Quest, 'n', ignore.case = TRUE) == 1){
        continue = TRUE
      }else if (regexpr(Quest, 'y', ignore.case = TRUE) == 1){
        st_write(joined, shape, delete_layer = TRUE)
      }
    }
  setwd(my.dir)
  }
}

#write.csv(processed_data, "processed_data.csv")
