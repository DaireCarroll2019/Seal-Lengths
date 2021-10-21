#Daire Carroll, Gothenburg University, 2021
#a function to measure the curved length of polygons representing seals

#########################################

if(!require("sf")){
  install.packages("sf",dependencies = TRUE)
}else if(!require("smoothr")){
  install.packages("smoothr",dependencies = TRUE)
}else if(!require("reshape2")){
  install.packages("reshape2",dependencies = TRUE)
}else if(!require("lwgeom")){
  install.packages("lwgeom",dependencies = TRUE)
}

#########################################

res = 100 #the resolution of the fit - used to densify the polygon for more points

clip_add = 30 #a small region at the end of each polygon is measured as a straight line to avoid overfitting (clip_add = number of coordinate points)

bu = 2 #the bandwidth value for smoothing - which removes or reduces the impact of seal limbs on measurments

#seals = st_transform(seals,23032) #format
seals = st_read("seal-pups_02-08-21_1241.shp")

processed_data[395,]

p = seals[39,]
curved_length(p,plt = TRUE)

#########################################

apex = function(poly){

	dist_pts = dist(cbind((st_coordinates(poly)[,1]),(st_coordinates(poly)[,2])),
	method = "euclidean", diag = TRUE, upper = TRUE
		)

	dist_pts = melt(as.matrix(dist_pts), varnames = c("row", "col"))

	loc = which(dist_pts$value == max(dist_pts$value)) 

	selection = data.frame(rbind(dist_pts[loc[1],],dist_pts[loc[length(loc)],]))

	pt1x = st_coordinates(poly)[,1][selection[1,1]]
	pt1y = st_coordinates(poly)[,2][selection[1,1]]
	pt1 = cbind(pt1x,pt1y)
	pt1 = st_point(pt1)

	pt2x = st_coordinates(poly)[,1][selection[1,2]]
	pt2y = st_coordinates(poly)[,2][selection[1,2]]
	pt2 = cbind(pt2x,pt2y)
	pt2 = st_point(pt2)

	pt1 = cbind(pt1x,pt1y)
	pt1 = st_point(pt1)
	pt2 = cbind(pt2x,pt2y)
	pt2 = st_point(pt2)

	pts = st_multipoint(rbind(pt1,pt2))

	return(pts)
	
}

midpoint = function(sf_lines = NULL){
	g = st_geometry(sf_lines)
	g_mids = lapply(g, function(x){
		coords = as.matrix(x)
		get_mids = function(coords){
			dist = sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
			dist_mid = sum(dist)/2
			dist_cum = c(0, cumsum(dist))
			end_index = which(dist_cum > dist_mid)[1]
			start_index = end_index - 1
			start = coords[start_index, ]
			end = coords[end_index, ]
			dist_remaining = dist_mid - dist_cum[start_index]
			mid = start + (end - start) * (dist_remaining/dist[start_index])
			return(mid)
		}
		return(get_mids(coords))
	})
	return(unlist(g_mids))
}

#########################################
#seven is two seals... consider adding watershed to overly large seals and see what they look like, otherwise exclude them

curved_length = function(p,plt){

	pol1 = p 
	pol2 = smooth(pol1, method = "ksmooth", smoothness = bu) #remove fins and flippers

	#plot(st_geometry(pol1), add = TRUE)
	plot(st_geometry(pol1))

	flipper_apex = apex(pol1)
	smooth_apex = apex(pol2)

  if(plt == TRUE){
	  plot(st_geometry(flipper_apex ), add = TRUE, col = "red")
  }
	
	buff = st_distance(flipper_apex,smooth_apex)

#	pol2 = st_buffer(pol2, buff)

	pline = st_cast(pol2, "LINESTRING")

	cuts = apex(pol2)

	parts = st_collection_extract(st_split(pline, cuts),"LINESTRING")

	lens = c()

	for(i in 1:length(row.names(parts))){
		lens[i] = length(st_coordinates(parts[i,]))	
	}

	if(length(lens)>2){
		equal_vals = which(lens == max(lens))
		if(length(equal_vals)>1){
			h1 = parts[equal_vals[1],]
			h2 = parts[equal_vals[2],] 
		}else{

			sections = c()

			for(i in 1:length(lens)){
				if(lens[i] == max(lens)){
					h1 = parts[i,]
				}else{
					sections[length(sections)+1] = i
				}
			}

			for(i in 1:(length(sections)-1)){
				h2 = st_union(parts[sections[i],],parts[sections[i+1],])
			} 	

		}

	}else{

		h1= parts[1,]
		h2 = parts[2,]

	}

	if(sum(st_coordinates(h1)[1,1:2]) != sum(st_coordinates(h2)[1,1:2])){
		XY = cbind(rev(st_coordinates(h1)[,1]),rev(st_coordinates(h1)[,2]))
		XY = st_linestring(x = XY, dim = "XY")
		XY = st_sfc(XY)
		h1 = st_set_geometry(h1,XY)
	}

	h1_l = length(st_coordinates(h1)[,1])
	h2_l = length(st_coordinates(h2)[,1])
	h1 = densify(h1,n = h2_l)
	h2 = densify(h2,n = h1_l)

	while(res*clip_add>(length(st_coordinates(h1)[,1])-res*clip_add)){
	  h1 = densify(h1,n = 2)
	  h2 = densify(h2,n = 2)
	}
	  
#	plot(st_geometry(h1))
#	plot(st_geometry(h2),add = TRUE, col = "green")
	if(plt == TRUE){
	  plot(st_geometry(h1),add = TRUE, col = "red")
	  plot(st_geometry(h2),add = TRUE, col = "green")
	}
	
	spine = data.frame(matrix(ncol = 2, nrow = 0))

	for(i in seq(res*clip_add,(length(st_coordinates(h1)[,1])-res*clip_add),by = res)){ #write this into the loop for checking for over estimation of length to geal with curve back on self bug

		tryCatch({
			pts = st_multipoint(rbind(st_coordinates(h1)[i,][1:2],st_coordinates(h2)[i,][1:2]))
			l2 = st_cast(pts, "MULTILINESTRING")
			l2 = st_sfc(l2)

			if(plt == TRUE){
			  plot(st_geometry(l2), add = TRUE, col = "red")
			}
			
			mp = midpoint(l2)
			
			if(plt == TRUE){
			  plot(st_geometry(st_point(mp)), add = TRUE, pch = "X")
			}
			
			spine = rbind(spine,  mp)
		}, error=function(e){cat("\n")})
	}

	dist(rbind(smooth_apex[1,]),rbind(spine[1,]),
	     method = "euclidean"
  )
	
	m = rbind(spine[length(spine[,1]),])
	m1 = smooth_apex[1,]
	m2 = smooth_apex[2,]
	d_pts1 = dist(rbind(m,m1), method= "euclidean")
	d_pts2 = dist(rbind(m,m2), method= "euclidean")
	
	if(d_pts2>d_pts1){
	  spine = rbind(spine,  smooth_apex[1,]) 
	  spine = rbind(smooth_apex[2,],  spine)
	}else{
	  spine = rbind(spine,  smooth_apex[2,])
	  spine = rbind(smooth_apex[1,],  spine)
	} #HERE
	
	colnames(spine) = c("X","Y")

	spine_line = data.matrix(spine)
	spine_line = st_linestring(spine_line)
	len = st_length(spine_line) + 2*buff 

	plot(st_geometry(spine_line), add = TRUE, col = "blue")
	
	st_crs(h1) = st_crs(pol1)
	st_crs(h2) = st_crs(pol1)
	
	spine_pts = st_as_sf(spine, coords = c("X","Y"))
	st_crs(spine_pts) = st_crs(pol1)
	
	spine_perimiter1 = st_distance(h1,st_as_sf(spine_pts), by_element = TRUE)
	spine_perimiter2 = st_distance(h2,st_as_sf(spine_pts), by_element = TRUE)
	
	if(max(spine_perimiter1) > max(spine_perimiter2) || max(spine_perimiter1) == max(spine_perimiter2)){
	  maxd = which(spine_perimiter1 == max(spine_perimiter1))
	  width = as.numeric(2*spine_perimiter1[maxd])
	  max_w = st_geometry(spine_pts)[maxd]

	  if(plt == TRUE){
      plot(max_w, add = TRUE, col = "black")
	  }
	  
	}else if(max(spine_perimiter1) < max(spine_perimiter2)){
	  maxd = which(spine_perimiter2 == max(spine_perimiter2))
	  width = as.numeric(2*spine_perimiter2[maxd])
	  max_w = st_geometry(spine_pts)[maxd]

	  if(plt == TRUE){
	    plot(max_w, add = TRUE, col = "black")
	  }

	}

	pline = st_cast(p,"LINESTRING")
	width2 = as.numeric(st_distance(max_w, pline)*2) #currently obsolete - the distance from the "waist" of the seal to the unsmoothed polygon, the issue of the flippers returns, functionally this code makes very little difference to the plots :) 
#  print(paste("width2",width2))
#	if(width2>width){
#    print("controling for smoothing")
#	  width = width2
#	}

  dims = c(len,width,pol1$lat,pol1$lon)
	
	return(dims)

}




