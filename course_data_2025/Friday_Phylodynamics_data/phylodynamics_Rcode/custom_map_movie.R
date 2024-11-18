# CUSTOM FUNCTIONS FOR MOVIE PLOTS - uses BEAST mcc trees - you must use read_MCC_tree.R
# S J Lycett
# this version 30 Jun 2022
# 3 July 2023 - option to exclude branches leading to certain tips (ee=excludedTips) in pts_and_hpds_at_time
# 14 Oct 2023 - included North pole settings
# 8 Nov 2023 - more fixed for ellipses

library(maps)
library(mapdata)
library(mapproj)

# needed to fit HPDs to ellipses
library(conicfit)

# 26 Nov 2021 - added some na.rm=TRUE
# 8 Sept 2022 - added another catch for ellipDirect
# 4 Jan 2023 - added another catch for ellipDirect
# 8 Nov 2023
fit_HPDs_to_standard <- function( tr=tr, npts=50, ltol=0.005 ) {
  # fit hpds to standard
  hpds <- vector("list",length(tr$lat80))
  for (j in 1:length(tr$lat80)) {
    #xc    <- mean(tr$lat80[[j]])
    #yc    <- mean(tr$lon80[[j]])
    ytemp  <- tr$lat80[[j]] #- xc #tr$latlon[j,1]
    xtemp  <- tr$lon80[[j]] #- yc #tr$latlon[j,2]
    
    xgood  <- (length(xtemp)>3)
    ygood  <- (length(ytemp)>3)
    if (xgood) xgood <- xgood & ( abs(max(xtemp,na.rm=TRUE)-min(xtemp,na.rm=TRUE))>ltol )
    if (ygood) ygood <- ygood & ( abs(max(ytemp,na.rm=TRUE)-min(ytemp,na.rm=TRUE))>ltol )
    
    if (xgood & ygood) {
      kk     <- which(is.finite(xtemp) & is.finite(ytemp))
      xy		 <- cbind(xtemp[kk],ytemp[kk])
      
      # 4 Jan 2023 & 8 Nov 2023
      nk     <- length(kk)
      if (nk>3) {
        # duplicated last entries cause issues with EllipseDirectFit
        if (all(xy[nk,]==xy[nk-1,]) | (all(xy[nk,]==xy[1,])) ) {
          xy <- xy[1:(nk-1),]
        }
	}
      ellipDirect  <- EllipseDirectFit(xy)
      if (all(is.na(ellipDirect))) {
        # small circle only
        ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),ltol,ltol,0)
      } else {
        ellipDirectG <- AtoG(ellipDirect)$ParG
      }
    } else {
      if (!xgood) {
        if (!ygood) {
          # small circle only
          ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),ltol,ltol,0)
        } else {
          # has y range only
          yrange <- (max(ytemp,na.rm=TRUE)-min(ytemp,na.rm=TRUE))/2
          ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),ltol,yrange,0)
        }
      } else {
        # has x range only
        xrange <- (max(xtemp,na.rm=TRUE)-min(xtemp,na.rm=TRUE))/2
        ellipDirectG<- c(mean(xtemp,na.rm=TRUE),mean(ytemp,na.rm=TRUE),xrange,ltol,0)
      }
    }
    
    
    xyFitted    <- calculateEllipse(ellipDirectG[1], ellipDirectG[2], 
                                    ellipDirectG[3], ellipDirectG[4], 
                                    180/pi*ellipDirectG[5], steps=npts)
    
    
    #plot(xyFitted[,1],xyFitted[,2], type="l", col="magenta")
    #points(mean(xtemp),mean(ytemp),pch=21,bg="red")
    #lines( c(mean(xtemp),mean(xtemp)), c(min(ytemp),max(ytemp)), col="blue")
    #lines( c(min(xtemp),max(xtemp)), c(mean(ytemp),mean(ytemp)), col="blue")
    #if (length(xtemp)==length(ytemp) & length(xtemp)>1) {
    #	polygon(xtemp,ytemp)			
    #} else {
    #	if (length(xtemp)==1) {
    #		points( array(xtemp, length(ytemp)), ytemp)
    #	}
    #	if (length(ytemp)==1) {
    #		points( xtemp, array(ytemp, length(xtemp)) )
    #	}
    #}
    
    hpds[[j]]	<- xyFitted
  }
  tr$hpds <- hpds
  
  # hpds differences
  hpd_diffs <- vector("list", length(tr$edge.length))
  for (j in 1:length(tr$edge.length)) {
    fromHPD <- hpds[[ tr$edge[j,1] ]]
    toHPD   <- hpds[[ tr$edge[j,2] ]]
    hpd_diffs[[j]] <- toHPD-fromHPD
  }
  tr$hpd_diffs <- hpd_diffs
  
  
  return( tr )
}


# 17 mar 2019
addFullTraitSet <- function(tr=tr, propIndex=1) {
  # see read_MCC_tree.R for getFullTraitSet function
  fullset 	<- getFullTraitSet(tr=tr, propIndex=propIndex)
  tr$fullset 	<- fullset
  return( tr )
}

# 17 mar 2019
interpolate_discrete_element <- function( fromSet=fromSet, toSet=toSet, fractTime=fractTime) {
  ff 	 <- matrix( rep(fractTime,length(fromSet[1,])), length(fractTime), length(fromSet[1,]) )
  midSet <- (fromSet*(1-ff)) + (toSet*ff)
  return( midSet )
}

# returns nodes and hpds suitable for plotting at a time point
# interpolates between branches which cut the time point
# optionally return the interpolated discrete trait if addFullTraitSet is used
# 3 July 2023 - option to exclude branches leading to certain tips (ee=excludedTips)
pts_and_hpds_at_time <- function(timePt, tr=tr, ii=c(), xlim=c(-180,180), ylim=c(-60,80), excludedTips=c()) {
  
  # latlon -2 = x, -1 = y
  fromY		<- tr$latlon[tr$edge[,1],1]
  fromX 	<- tr$latlon[tr$edge[,1],2]
  toY		<- tr$latlon[tr$edge[,2],1]
  toX 		<- tr$latlon[tr$edge[,2],2]
  fromTime 	<- tr$nodeTimes[tr$edge[,1]]
  toTime   	<- tr$nodeTimes[tr$edge[,2]]
  branchTime  <- (toTime-fromTime)
  branchX	<- (toX-fromX)
  branchY	<- (toY-fromY)
  
  # 16 mar 2019
  hpds 		<- tr$hpds
  hpd_diffs 	<- tr$hpd_diffs
  
  # identify branches
  if (length(ii)==0) {
    ii 		<- which( toTime >= timePt & fromTime <= timePt  &
                     toX >= xlim[1] & toX <= xlim[2] & toY >= ylim[1] & toY <= ylim[2] )
  }
  if (length(excludedTips)>0) {
	ex_tips <- match(excludedTips,tr$tip.label)
	ex_tips <- ex_tips[which(is.finite(ex_tips))]
	if (length(ex_tips)>0) {
		ee <- match(ex_tips,tr$edge[,2])
		print(paste("Excluding",length(ee),"branches to tips"))
		ii <- setdiff(ii,ee)
	}
  }
  
  fractTime 	<- 1-(toTime-timePt)/(branchTime)
  
  if (length(ii)>0) {
    if (length(ii)==1) {
      ii <- c(ii,ii)
    }
    # get mcc pts
    x_pos		<- (fractTime*branchX + fromX)[ii]
    y_pos		<- (fractTime*branchY + fromY)[ii]
    
    pts 		<- list(x=x_pos,y=y_pos)
    
    # hpds
    pts_polys <- vector("list",length(ii))
    for (j in 1:length(ii)) {
      fromNode	  <- tr$edge[ii[j],1]
      temp_hpd_diff <- hpd_diffs[[ii[j]]]
      temp_hpd	  <- fractTime[ii[j]]*temp_hpd_diff + hpds[[fromNode]]
      pts_polys[[j]] <- temp_hpd
    }
  } else {
    pts       <- NULL
    pts_polys <- NULL
  }
  
  if ( any(attributes(tr)$names=="fullset") ) {
    # interpolate the selected discrete trait
    # use for colouring
    if (length(ii)>0) {
      fromSet <- tr$fullset[tr$edge[ii,1],]
      toSet	  <- tr$fullset[tr$edge[ii,2],]
      midSet  <- interpolate_discrete_element(fromSet=fromSet, toSet=toSet, fractTime=fractTime[ii])
      interpol_trait <- colnames(midSet)[apply(midSet, 1, which.max)]
    
      return( list(pts=pts, pts_polys=pts_polys, 
                 ii=ii, fractTime=fractTime[ii], 
                 interpol_trait=interpol_trait, midSet=midSet) )
    } else {
      return( list(pts=pts, pts_polys=pts_polys, ii=ii, fractTime=fractTime[ii]) )  
    }
  } else {
    return( list(pts=pts, pts_polys=pts_polys, ii=ii, fractTime=fractTime[ii]) )
  }
}

# 17 mar 2019; 5 Nov 2021; 8 Sept 2022; 15 Oct 2023
plot_at_time <- function(tpt=tpt, timePt=timePt, tr=tr, xlim=c(-180,180), ylim=c(-60,80), 
                         show.hpds=TRUE, solid.pts=TRUE, show.legend=TRUE, new.plot=TRUE,
                         timelegpos="bottomright",timeFormat="ddmmyy",
                         legpos="bottomleft",lcex=1,
                         fcol="white", bgcol="grey90", bdcol="grey70", fill=TRUE,
                         use.fullset.cols=TRUE,
                         h1=0, h2=0,
                         s1=0.3, s2=0.7,
                         b1=0.9, b2=0.7,
                         t1=0.25, t2=0.75,
                         ppch=21,
                         useWorldHires=TRUE,
                         utraits=colnames(tr$fullset),
                         tcols=get_BEAST_cols( length(utraits), sat=s1, bright=b1, transparency=t1 ),
                         tcols2=get_BEAST_cols( length(utraits), sat=s2, bright=b2, transparency=t2 )
                         ) {
  
  if ( length(tpt)== 0 ) {
    tpt 	<- pts_and_hpds_at_time(timePt, tr=tr, xlim=xlim, ylim=ylim)
  }
  
  npcol <- length(tpt$pts_polys)
  if (length(npcol) < 1) {
    npcol <- 1
  }
  pcol 	<- array(hsv(h1, s1, b1, t1), npcol)
  pcol2 <- array(hsv(h2, s2, b2, t2), npcol)
  
  if (new.plot) {
    if (useWorldHires) {
      map("worldHires", xlim=xlim, ylim=ylim, 
          col=fcol, fill=fill, border=bdcol, bg=bgcol)
    } else {
      map("world", xlim=xlim, ylim=ylim, 
          col=fcol, fill=fill, border=bdcol, bg=bgcol)
    }
  }
  
  if ( use.fullset.cols & any(attributes(tr)$names=="fullset") ) {
    #utraits <- colnames(tr$fullset)
    #tcols   <- get_BEAST_cols( length(utraits), sat=s1, bright=b1, transparency=t1 )
    #tcols2  <- get_BEAST_cols( length(utraits), sat=s2, bright=b2, transparency=t2 )
    
    tinds	  <- match(tpt$interpol_trait, utraits)
    if (length(tinds)==0) {
      tinds <- length(utraits)
    }
    pcol	  <- tcols[tinds]
    pcol2	  <- tcols2[tinds]
    
    if (show.legend & legpos!="x") {
      legend(legpos, paste(utraits), pch=ppch, pt.bg=tcols2, bty="n", lcex=lcex)
    }
  }
  
  if (show.hpds & length(tpt$pts_polys)>0) {
    for (j in 1:length(tpt$pts_polys)) {
      polygon(tpt$pts_polys[[j]], col=pcol[j], bg=pcol[j], border=FALSE)
    }
  }
  
  if (length(tpt$pts)>0) {
    if (solid.pts) {
      if (ppch<21) {
        points(tpt$pts, pch=ppch, col=pcol2)
      } else {
        points(tpt$pts, pch=ppch, bg=pcol2)
      }
    } else {
      points(tpt$pts, col=pcol2)
    }
  }
  
  if (show.legend & timelegpos!="x") {
    if (timeFormat=="ddmmyy" | timeFormat=="mmyy") {
      timeTxt <- invertDecimalDate(timePt,ddmmyy=TRUE)
      if (timeFormat=="mmyy") {
        timeTxt <- substring(timeTxt, 4)
      }
      legend(timelegpos, timeTxt, bty="n", pch=NA)
    } else {
      legend(timelegpos, format(timePt,digits=7), bty="n", pch=NA)
    }
  }
}


plot_mcc_tree_with_hpds <- function(tr=tr, xlim=c(-180,180), ylim=c(-60,80), 
                                    tlim=c(NA,NA),
                                    propIndex=0,
                                    show.hpds=TRUE, solid.pts=TRUE, show.legend=TRUE, 
                                    timelegpos="bottomright",
                                    legpos="bottomleft",
                                    useWorldHires=FALSE,new.plot=TRUE,
                                    fcol="white", bgcol="grey90", bdcol="grey70", fill=TRUE,
						                        h1=0, h2=0.6,
                                    s1=0.3, s2=0.7,
                                    b1=0.9, b2=0.7,
                                    t1=0.25, t2=0.75,
						                        utraits=tr$uprops[[propIndex]],
						                        tcols=c(-1), tcols2=c(-1),lcex=1) {
  
  pcol  <- array(hsv(h1, s1, b1, t1), length(tr$hpds))
  pcol2 <- array(hsv(h2, s2, b2, t2), length(tr$hpds))
  ecol2 <- array(hsv(0, 0, 0.2, t2), length(tr$edge))
  
  if (new.plot) {
    if (useWorldHires) {
      map("worldHires", xlim=xlim, ylim=ylim, 
        col=fcol, fill=fill, border=bdcol, bg=bgcol)
    } else {
      map("world", xlim=xlim, ylim=ylim, 
          col=fcol, fill=fill, border=bdcol, bg=bgcol)
    }
  }
  
  if (propIndex>0) {
    #utraits <- tr$uprops[[propIndex]]
    if (tcols[1]==-1 & tcols2[1]==-1) {
     	tcols   <- get_BEAST_cols( length(utraits), sat=s1, bright=b1, transparency=t1 )
     	tcols2  <- get_BEAST_cols( length(utraits), sat=s2, bright=b2, transparency=t2 )
    }
    
    tinds	  <- match(tr$props[,propIndex], utraits)
    pcol	  <- tcols[tinds]
    pcol2	  <- tcols2[tinds]
    einds	  <- match(tr$props[tr$edge[,1],propIndex], utraits)
    ecol2	  <- tcols2[einds]
    
    if (show.legend) {
      legend(legpos, paste(utraits), pch=21, pt.bg=tcols2, bty="n", cex=lcex)
    }
  }

  #8 july 2022 - temporal cols
  if (propIndex==-1) {
	if (is.na(tlim[1])) {
		tt1 <- floor(min(tr$nodeTimes))
		tt2 <- ceiling(max(tr$nodeTimes))
	} else {
		tt1 <- tlim[1]
		tt2 <- tlim[2]
	}
	cpts <- 100
	timecols  <- hsv( h2-((h2-h1)*(0:(cpts-1))/cpts), s1, b1, t1)
	timecols2 <- hsv( h2-((h2-h1)*(0:(cpts-1))/cpts), s2, b2, t2)
	normT     <- ceiling(cpts*(tr$nodeTimes-tt1)/(tt2-tt1))
	normT[which(normT<1)]    <- 1
	normT[which(normT>cpts)] <- cpts
	pcol <- timecols[normT]
	pcol2<- timecols2[normT]
	ecol2<- timecols2[normT[tr$edge[,2]]]

	ll    <- c(1,25,50,75,100)
	lcols <- timecols2[ll]
	lvals <- c(tt1,tt1+(tt2-tt1)*0.25,tt1+0.5*(tt2-tt1),tt1+0.75*(tt2-tt1),tt2)
	if (show.legend) {
		legend(legpos,format(lvals,digits=3),lty=1,col=lcols,bty="n",cex=lcex)
	}
  }
  
  tr$ntips <- length(tr$tip.label)
  tinds    <- 1:tr$ntips
  ninds    <- (tr$ntips+1):length(tr$latlon[,1])
  fromY		<- tr$latlon[tr$edge[,1],1]
  fromX 	<- tr$latlon[tr$edge[,1],2]
  toY		<- tr$latlon[tr$edge[,2],1]
  toX 		<- tr$latlon[tr$edge[,2],2]
  
  if (!is.na(tlim[1])) {
    ok_nodes <- which(tr$nodeTimes >= tlim[1] & tr$nodeTimes <= tlim[2])
    tinds	   <- intersect(tinds, ok_nodes)
    ninds	   <- intersect(ninds, ok_nodes)
    ok_edges <- which(tr$nodeTimes[tr$edge[,2]] >= tlim[1] & 
                        tr$nodeTimes[tr$edge[,1]] <= tlim[2] )
  } else {
    ok_nodes <- c(tinds,ninds)
    ok_edges <- 1:length(tr$edge[,1])
  }
  
  #for (i in 1:length(tr$hpds)) {
  if (show.hpds) {
    for (j in length(ok_nodes):1) {
      i <- ok_nodes[j]
      polygon(tr$hpds[[i]], col=pcol[i], bg=pcol[i], border=FALSE)
    }
  }
  
  arrows( 	fromX[ok_edges], 	fromY[ok_edges], 
           toX[ok_edges], 	toY[ok_edges], 
           length=0.1, angle=15, col=ecol2[ok_edges])
  
  points(tr$latlon[tinds,2], tr$latlon[tinds,1], pch=21, bg=pcol2[tinds])
  points(tr$latlon[ninds,2], tr$latlon[ninds,1], pch=21, col=pcol2[ninds])
  
}

##############################################################
# North pole settings
# 14 Oct 2023

# 14 oct 2023 - sub function to reverse order of latlon coordinates used by switch_latlon
# useful if these are really x,y (northpole) and you want to use the above functions
switch_list_el <- function(el) {
  return( el[,c(2,1)] )
}

# 14 oct 2023 - function to reverse order of latlon coordinates
# useful if these are really x,y (northpole) and you want to use the above functions
# e.g. north_pole_settings(); plot_mcc_tree_with_hpds(tr=switch_latlon(tr),new.plot=FALSE)
switch_latlon <- function(tr) {
  orig_latlon <- tr$latlon
  orig_lat80  <- tr$lat80
  orig_lon80  <- tr$lon80
  orig_hpds   <- tr$hpds
  orig_hpds_diffs <- tr$hpd_diffs
  
  new_latlon     <- switch_list_el(orig_latlon)
  new_lat80      <- orig_lon80
  new_lon80      <- orig_lat80
  new_hpds       <- lapply(orig_hpds,switch_list_el)
  new_hpds_diffs <- lapply(orig_hpds_diffs,switch_list_el)
  
  new_tr <- tr
  new_tr$latlon <- new_latlon
  new_tr$lat80  <- new_lat80
  new_tr$lon80  <- new_lon80
  new_tr$hpds   <- new_hpds
  new_tr$hpd_diffs <- new_hpds_diffs
  
  return(new_tr)
}

north_pole_settings <- function(xlim=c(-180,180),ylim=c(-63,90),
                                ringlat=c(60,40,20,0,-20,-40,-63),
                                fcol="grey70",bgcol="white",bdcol=NA,fill=TRUE,
                                az_col="grey50",last_az_col=az_col,
                                np_col="grey30", show_map=TRUE) {
  
  projection 	<- "azequalarea" # "mercator" #"azequalarea"
  orientation <- c(90,0,0) 
  #xlim        <- c(-180,180)
  #ylim        <- c(-63,90)
  azs         <- seq(xlim[1],xlim[2],1)
  nring       <- length(ringlat)
  rings       <- matrix(0,nring,length(azs))
  cs          <- vector("list",nring)
  for (j in 1:nring) {
    rings[j,] <- ringlat[j]
    cs[[j]] <- mapproject(azs,rings[j,],projection=projection, orientation=orientation)
  }
  
  #fcol=    "grey70" #"white"
  #bgcol=   "white"  # grey80" #"grey90"
  #bdcol=   NA#"grey70"
  #fill=    TRUE
  #az_col<- "grey50"#"grey90"
  #np_col<- "grey30"
 
  # show map
  if(show_map) {
    #op <- par(mar=c(1,1,1,1))
    map("world",projection=projection, orientation=orientation,
        xlim=xlim,ylim=ylim,
        col=fcol, fill=fill, border=bdcol, bg=bgcol)
    points(0,0,pch=3,col=np_col)
    for (j in 1:nring) {
      if (j==nring) {
        lines(cs[[j]],col=last_az_col)
      } else {
        lines(cs[[j]],col=az_col)
      }
    }
    #par(op) 
  }
  
}

