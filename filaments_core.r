
# Filament & Tendril analysis functions
# Mehmet Alpaslan

## Notes: For the time being working on G12 out to redshift z = 0.2. Main sample has no magnitude limit.
## G12cutxyz - Cartesian positions for ALL GALAXIES in G12 out to z = 0.2
## G12sampxyz - Cartesian positions for ALL GROUPS in G12 out to z = 0.2
## Link between these and cats using CATAID and GroupID.

# Some really quick commands.
# plot a whole cone in 3d

plotcone <- function(){plot3d(G12sampxyz$X,G12sampxyz$Y,G12sampxyz$Z,aspect=c(0.5,2,0.3),size=2,xlab="X",ylab="Y",zlab="Z")}
sourceAll <- function(){
	source("~/Dropbox/filaments/code/filaments_core.r")
	source("~/Dropbox/filaments/code/filaments_functions.r")
	source("~/Dropbox/filaments/tendrils/tendril_functions.r")
	source("~/Dropbox/filaments/code/walker.R")
	source("~/Dropbox/filaments/galprops/galprops_functions.r")
	source("~/Dropbox/filaments/alongFilament/alongFilament_functions.r")
	cat("Code loaded!","\n")
}

library(Cairo)
library(astro)
library(nnclust)
library(RColorBrewer)
library(fields)
library(hyper.fit)
source("~/Dropbox/filaments/code/walker.R")

##################################################################################

# sph2car (by ASGR) - sphericals to cartesians
sph2car <- function(long,lat,radius=1,deg=T){
	if (is.matrix(long) || is.data.frame(long)){
	if(ncol(long)==1){lat=lat}else{
	if(ncol(long)==2){lat=long[,2];long=long[,1]}else{
	if(ncol(long)==3){radius=long[,3];lat=long[,2];long=long[,1]}
	}
	}
	}
	if(deg){long=long*pi/180;lat=lat*pi/180}
	return=cbind(radius*cos(long)*cos(lat),radius*sin(long)*cos(lat),radius*sin(lat))
}

# car2sph (by ASGR) - cartesian to sphericals    
car2sph <- function(xyz){
       long=atan2(xyz[,2],xyz[,1])
       if(long>0){lat=atan2(xyz[,3],xyz[,2]/sin(long))}else{
	lat=atan2(xyz[,3],xyz[,1]/cos(long))}
       radius=sqrt(xyz[,1]^2+xyz[,2]^2+xyz[,3]^2)
       return=cbind(long,lat,radius)
}

# CosDist (by ASGR) - cosmological distance calculator
CosDist <- function(z,H0=71.0,OmegaM=0.27,OmegaL=0.73,OmegaK=0,age=F){
    temp=function(z,H0=71.0,OmegaM=0.27,OmegaL=0.73,OmegaK=0){
    Einv=function(z,OmegaM=0.27,OmegaL=0.73,OmegaK=0){1/sqrt(OmegaM*(1+z)^3+OmegaK*(1+z)^2+OmegaL)}
    CoDist=(299792.458/H0)*integrate(Einv,0,z,OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000)$value
    CoVol=(4/3)*pi*CoDist^3
    LumDist=(1+z)*CoDist
    AngDist=CoDist/(1+z)
    AngArcSec=AngDist*(pi/(180*60*60))
    if(age){
    Einvz=function(z,OmegaM=0.27,OmegaL=0.73,OmegaK=0){1/(sqrt(OmegaM*(1+z)^3+OmegaK*(1+z)^2+OmegaL)*(1+z))}
    HT=3.08568025e19/(H0*31556926)
    UniAge=HT*integrate(Einvz,0,Inf,OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000)$value
    zAge=HT*integrate(Einvz,0,z,OmegaM=OmegaM,OmegaL=OmegaL,OmegaK=OmegaK,subdivisions=1000)$value
    }
    if(age){
    return=c(z=z,CoDist=CoDist,LumDist=LumDist,AngDist=AngDist,AngArcSec=AngArcSec,CoVolGpc3=CoVol/1e9,HubTime=HT,UniAge=UniAge,TravelTime=zAge)
    }else{
    return=c(z=z,CoDist=CoDist,LumDist=LumDist,AngDist=AngDist,AngArcSec=AngArcSec,CoVolGpc3=CoVol/1e9)
    }
    }
    return=t(Vectorize(temp)(z,H0,OmegaM,OmegaL,OmegaK))
}

#KECorr (by ASGR) - K correction (first output)

KEcorr=function(z){
k=sum(kcorrvals*(z-0.2)^(0:4))
return=c(k,k-1.75*z)}

kcorrvals=c(0.20848,1.0226,0.52366,3.5902,2.3843)

#segvol by ASGR. Calculates sphere segment volume/area.

segvol=function(lorad,hirad,lodec,hidec,lora,hira){
lodec=lodec*pi/180
hidec=hidec*pi/180
h=(sin(hidec)-sin(lodec))
area=2*pi*((hira-lora)/360)*h
vol=(((4/3)*pi*hirad^3)*area/(4*pi))-(((4/3)*pi*lorad^3)*area/(4*pi))
return=list(vol=vol,area=area)
}

# overlay grid by MA. Overlays grid on a coneplot. Supply angle/distance limits & cosmology.

overlaygrid <- function(ralim,declim,zlim,H0=100,OmegaM=0.25,OmegaL=0.75, text=c(1,2,3),textsize = 0.5, textangle=c(0,0),radeclabels = F, rotangle = NA){

	if(is.na(rotangle)==T){rotangle = ((max(ralim)+min(ralim))/2) - 90}
	rasteps = seq(min(ralim),max(ralim),by=2)
	zsteps = seq(min(zlim),max(zlim),by=0.01)
	
	toplines = rotate3d((sph2car(rasteps,0,comovdist.los(z=max(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	bottomlines = rotate3d((sph2car(rasteps,0,comovdist.los(z=min(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	topright = rotate3d((sph2car(min(ralim),0,comovdist.los(z=max(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	topleft = rotate3d((sph2car(max(ralim),0,comovdist.los(z=max(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	bottomright = rotate3d((sph2car(min(ralim),0,comovdist.los(z=min(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	bottomleft = rotate3d((sph2car(max(ralim),0,comovdist.los(z=min(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	
	for(i in 1:length(rasteps)){
	lines(x=c(toplines[i,1],bottomlines[i,1]),y=c(toplines[i,2],bottomlines[i,2]),col=hsv(h=0,s=0,v=0.75,alpha=0.75),lwd=0.75)
	if(any(text==1 & i %% 2)){text(toplines[i,1],toplines[i,2]+3.75,rasteps[length(rasteps)+1-i],cex=textsize, srt=textangle[1])}
	}
	lines(x=c(topright[1],bottomright[1]),y=c(topright[2],bottomright[2]))
	lines(x=c(topleft[1],bottomleft[1]),y=c(topleft[2],bottomleft[2]))
	
	toparc = rotate3d((sph2car(seq(min(ralim),max(ralim),0.1),0,comovdist.los(z=max(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	bottomarc = rotate3d((sph2car(seq(min(ralim),max(ralim),0.1),0,comovdist.los(z=min(zlim),H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)
	
	for(i in 1:length(zsteps)){
	arc = rotate3d((sph2car(seq(min(ralim),max(ralim),0.1),0,comovdist.los(z=zsteps[i],H=H0,M=OmegaM))),(rotangle*pi/180),0,0,1)	
	lines(x=arc[,1],y=arc[,2],col=hsv(h=0,s=0,v=0.75,alpha=0.75),lwd=0.75)
	if(any(text==2 & i%%2 == 0)){text(head(arc[,1],1)+7.15,tail(arc[,2],1),zsteps[i],cex=textsize, srt=textangle[2])}
	}
	
	lines(x=toparc[,1],y=toparc[,2])
	lines(x=bottomarc[,1],y=bottomarc[,2])
	
	if(any(text==3)){text(x=(par('usr')[1]+par('usr')[2])/2,y=par('usr')[3]+10,label=paste(min(zlim), " < z < ", max(zlim), "\n",min(declim), " < Dec < ", max(declim),sep=""),cex=textsize, srt=textangle[2])}

	if(radeclabels==T){
		text(x=(par('usr')[1]+par('usr')[2])/2, y=par('usr')[4]-25,label="RA [deg]",cex=1.25)
		text(x=par('usr')[2]-45, y=(par('usr')[3]+par('usr')[4])/2-20, label="z",cex=1.25)
	}

}

##################################################################################

# groupit by ASGR; edited by MA. Takes input as a 2 column table containing from and to links, and spits out groups with group numbers. This only considers something a group if links go both ways -- bear this in mind.

# Comments by MA.

groupit <- function(ind){

	# Start by using table to generate a list of how many times the from column is repeated, and sort in descending order.

	table=table(ind[, 1])
	table=cbind(as.numeric(names(table)),as.numeric(table))
	table=table[order(table[,2],decreasing=T),]
	
	# Check if there are any pairs. If not, skip making grefs2.
	
	dopairs = length(ind[ind[,1] %in% table[table[,2]==1,1] & ind[,2] %in% table[table[,2]==1,1],])
	
	# grefs2 is designed to identify all pairs within the data. Do this by extracting everything in ind that appears once in table.
	if(dopairs > 0){
	
		grefs2=ind[ind[,1] %in% table[table[,2]==1,1] & ind[,2] %in% table[table[,2]==1,1],]
		grefs2=as.numeric(grefs2[rep(1:length(grefs2[,1]),each=2)+c(0,length(grefs2[,1]))])
		grefs2=unique(grefs2)
		if(length(grefs2)>0){
		grefs2=cbind(grefs2,rep(1:(length(grefs2)/2),each=2))
		}
    
	# Remove pairs from table.

		table=table[! table[,1] %in% grefs2[,1],]
	}else{grefs2=c(0,0)}
	
	# If data is all pairs, skip straight to the output.

	if(dim(table)[1] > 0){

		# Set up grefs and begin generating groups.
		
	    grefs = cbind(table[,1], 0)
	    print("Generating Groups")
	    gnum=1
	    while (any(grefs[, 2] == 0)) {
	        if (gnum%%1000 == 0) {
	            print(paste("Group", gnum, "with", length(which(grefs[, 
	                2] == 0)), "galaxies left"))
	        }
	        
	    # Generate groups from unassigned objects in grefs.
	        
	        tempref = min(which(grefs[, 2] == 0))
	        tempref = grefs[tempref, 1]
	        extract = ind[ind[, 1] == tempref, 2]
	        
	    # Define a group as the object from FROM and all objects in TO that it leads to.
	    
	        group = c(tempref, extract)
	        while (length(extract) > 0) {
	            extract = unique(ind[ind[, 1] %in% group & !ind[,2] %in% group, 2])
	            group = c(group, extract)
	        }
	        
	        grefs[grefs[, 1] %in% group, 2] = gnum
	        
	    # Once a group has been identified, remove it from ind and continue.
	    
	        ind = ind[!(ind[, 1] %in% group | ind[, 2] %in% group), ]
	        gnum = gnum + 1
	    }
	    if(length(grefs2)>0 & dopairs > 0){
	   grefs2[,2]=grefs2[,2]+max(grefs[,2])
	   grefs=rbind(grefs,grefs2)
	    }

	}

	if(dim(table)[1]==0){grefs = grefs2}

   colnames(grefs)=c('groupID','filID')
   return=grefs
    }


#_________________________________________#
#_________________________________________#

# createSample
# This function is designed to create a volume limited sample of groups for a maximum redshift. It does this by calculating the faintest possible absolute magnitude that can be observed at that redshift. It removes any galaxies in G3CGal that are fainter. Then it removes them from the group catalogue; and any groups that have fewer than 2 members left are removed. Outputs are galaxies and groups. Inputs are just redshift.

createSample <- function(z, groupcat, galcat, appmaglim = 19.4, plot=F){
	
# Grab k-correction value for the redshift, then the distance-modulus.

	k = KEcorr(z)[1]
	DM = 5*log10(CosDist(z,H0=100,OmegaM=0.25,OmegaL=0.75)[,'LumDist']*1E5)
	M = appmaglim - DM - k

	print(M)

# Keep any galaxies fainter than M.

	kfunction = approxfun(x=seq(0,0.6,by=0.0005),sapply(seq(0,0.6,by=0.0005),KEcorr)[1,])
	gal_k = sapply(galcat[,4],kfunction)
	
	GalAbsMag = galcat$Rpetro - galcat$DM_100_25_75 - gal_k

	keeprows = which(GalAbsMag <= M)

	if(colnames(groupcat)[1] == "GroupID"){galid = grep("GroupID",colnames(galcat));groupid=grep("GroupID",colnames(groupcat))}else{galid = grep("HaloID",colnames(galcat));groupid = grep("HaloID",colnames(groupcat))}

# Loop through all groups and remove galaxies fainter than M. If the group has fewer than 2 members left, discard group.

	discardcat = {}

	for(i in 1:length(groupcat[,1])){

		groupgals = cbind(galcat[which(galcat[,galid]==groupcat[i,groupid]),],GalAbsMag[which(galcat[,galid]==groupcat[i,groupid])])
		colnames(groupgals) = c(colnames(galcat),"AbsMag")

		keepgals = which(groupgals[,'AbsMag'] <= M)

# Check how many groups are kept; if less than 2, discard group. Discardgroup 1 means discarded, 0 means kept.

		if(length(keepgals) < 2){discardgroup = 1}else{discardgroup = 0}

		discardcat = rbind(discardcat,c(groupcat[i,groupid],discardgroup))

	}

	if(plot){

		magplot(groupcat$Zfof[which(discardcat[,2]==0)],groupcat$TotFluxProxy[which(discardcat[,2]==0)],pch='.',log='y',xlab='z',ylab=expression(paste("Group luminosity (L"["ʘ"]," h"^{"-2"},")",sep="")),xlim=c(0,0.5),ylim=c(5E9,1E13))
		abline(v = z, col='red')
		shade(x=c(z,par('usr')[2]),y1=c(par('usr')[3],par('usr')[3]),y2=c(10^par('usr')[4],10^par('usr')[4]),col=hsv(alpha=0.35))
		label('topleft',paste("Kept = ", length(which(discardcat[,2]==0)),'\n',"Discarded = ", length(which(discardcat[,2]==1)),sep=''))
		label('bottomleft',length(which(groupcat$Zfof <= z & discardcat[,2]==0)))
		label('bottomright',length(which(groupcat$Zfof > z & discardcat[,2]==0)))

	}

	return(discardcat)

}

#_________________________________________#
#_________________________________________#

# limitCone
# This is a simple function whose purpose is to remove all groups within x Mpc from the survey edge and spit out a modified catalogue. Catalogue should just be the group catalogue., cone limits (G09, G12, etc.) and the cutting distance. To calculate planes, use one point as origin, and other two as edges of plane at low redshift. 

limitCone <- function(groupcat,ralim,declim,zlim,cut){

	# Start by identifying which cone we're dealing with. We want to rotate everything so the cone originates at (0,0,0) and faces upwards, so we need a rotation angle, which is determined by the cone.

	if(min(ralim) > 120 & max(ralim) < 145){rotangle=45;field="G09"}
	if(min(ralim) > 170 & max(ralim) < 200){rotangle=90;field="G12"}
	if(min(ralim) > 210 & max(ralim) < 250){rotangle=127.5;field="G15"}

	# Make a 4 column cartesian coordinate data frame. Convert RA, Dec, z to comoving cartesians & rotate.

	zcol = grep("IterCenZ",colnames(groupcat))+1

	rows = which(groupcat$IterCenRA >= min(ralim) & groupcat$IterCenRA <= max(ralim))
	temp = sph2car(groupcat$IterCenRA[rows],groupcat$IterCenDec[rows],CosDist(groupcat[rows,zcol],100,0.25,0.75)[,'CoDist'])
	temp = rotate3d(temp,(rotangle*pi/180),0,0,1)
	temp = cbind(temp,groupcat[rows,1])

	cat = temp;colnames(cat) = c("X","Y","Z","GroupID");cat = as.data.frame(cat)

	# Now proceed.

	right = rotate3d((sph2car(min(ralim),declim,comovdist.los(z=as.numeric(quantile(zlim,c(0.25))),H=100,M=0.25))),(rotangle*pi/180),0,0,1)
	left = rotate3d((sph2car(max(ralim),declim,comovdist.los(z=as.numeric(quantile(zlim,c(0.25))),H=100,M=0.25))),(rotangle*pi/180),0,0,1)
	bottom = rotate3d((sph2car(ralim,min(declim),comovdist.los(z=as.numeric(quantile(zlim,c(0.25))),H=100,M=0.25))),(rotangle*pi/180),0,0,1)
	top = rotate3d((sph2car(ralim,max(declim),comovdist.los(z=as.numeric(quantile(zlim,c(0.25))),H=100,M=0.25))),(rotangle*pi/180),0,0,1)

	# Now define planes describing left, top, bottom, right of cone. For each point above, calculate the normal vector.

	# Take cross products.

	topn = c((top[1,2]*top[2,3])-(top[1,3]*top[2,2]),(top[1,3]*top[2,1])-(top[1,1]*top[2,3]),(top[1,1]*top[2,2])-(top[1,2]*top[2,1]))
	bottomn = c((bottom[1,2]*bottom[2,3])-(bottom[1,3]*bottom[2,2]),(bottom[1,3]*bottom[2,1])-(bottom[1,1]*bottom[2,3]),(bottom[1,1]*bottom[2,2])-(bottom[1,2]*bottom[2,1]))
	leftn = c((left[1,2]*left[2,3])-(left[1,3]*left[2,2]),(left[1,3]*left[2,1])-(left[1,1]*left[2,3]),(left[1,1]*left[2,2])-(left[1,2]*left[2,1]))
	rightn = c((right[1,2]*right[2,3])-(right[1,3]*right[2,2]),(right[1,3]*right[2,1])-(right[1,1]*right[2,3]),(right[1,1]*right[2,2])-(right[1,2]*right[2,1]))

	# Calculate d.

	topd = -(topn[1]*top[1,1]) - (topn[2]*top[1,2]) - (topn[3]*top[1,3])
	bottomd = -(bottomn[1]*bottom[1,1]) - (bottomn[2]*bottom[1,2]) - (bottomn[3]*bottom[1,3])
	leftd = -(leftn[1]*left[1,1]) - (leftn[2]*left[1,2]) - (leftn[3]*left[1,3])
	rightd = -(rightn[1]*right[1,1]) - (rightn[2]*right[1,2]) - (rightn[3]*right[1,3]) 

	# Now work out distance from each plane for all groups. Compare minimum distance to the cut value and choose to reject/keep.

	cutcat = {}
	for(i in 1:length(cat[,1])){

		dist = c(abs((topn[1]*cat[i,1]) + (topn[2]*cat[i,2]) + (topn[3]*cat[i,3]) + topd)/(sqrt(topn[1]^2+topn[2]^2+topn[3]^2)), abs((bottomn[1]*cat[i,1]) + (bottomn[2]*cat[i,2]) + (bottomn[3]*cat[i,3]) + bottomd)/(sqrt(bottomn[1]^2+bottomn[2]^2+bottomn[3]^2)), abs((leftn[1]*cat[i,1]) + (leftn[2]*cat[i,2]) + (leftn[3]*cat[i,3]) + leftd)/(sqrt(leftn[1]^2+leftn[2]^2+leftn[3]^2)), abs((rightn[1]*cat[i,1]) + (rightn[2]*cat[i,2]) + (rightn[3]*cat[i,3]) + rightd)/(sqrt(rightn[1]^2+rightn[2]^2+rightn[3]^2)))
		if(min(dist) >= cut){cutcat=rbind(cutcat,cat[i,])}

	}

	# Now do a radial cut.

	zcutcat = groupcat[groupcat[,1] %in% cutcat[,4],zcol]
	distcutcat = comovdist.los(z=zcutcat,H=100,M=0.25)
	limit = comovdist.los(z=max(zlim),H=100,M=0.25) - cut
	cutcat = cutcat[which(distcutcat <= limit),]

	# Want to add a special marker to all galaxies that are within the cut region. So add a new column onto the cat variable, where 0 means included, 1 means cut.

	temp1 = cbind(cat[(which(!cat[,1] %in% cutcat[,1])),],1);colnames(temp1) = c("X","Y","Z","GroupID","Cut")
	temp2 = cbind(cat[(which(cat[,1] %in% cutcat[,1])),],0);colnames(temp2) = c("X","Y","Z","GroupID","Cut")
	cutcat = rbind(temp1,temp2)

	# Reorder catalogue to go by GroupID.

	cutcat = cutcat[order(cutcat[,1]),]

	# plotcone()
	# points3d(top,size=5,col='blue')
	# points3d(bottom,size=5,col='green')
	# points3d(right,size=5,col='red')
	# points3d(left,size=5,col='yellow')
	# planes3d(a=topn[1],b=topn[2],c=topn[3],d=topd,col='blue',alpha=0.25)
	# planes3d(a=bottomn[1],b=bottomn[2],c=bottomn[3],d=bottomd,col='green',alpha=0.25)
	# planes3d(a=leftn[1],b=leftn[2],c=leftn[3],d=leftd,col='yellow',alpha=0.25)
	# planes3d(a=rightn[1],b=rightn[2],c=rightn[3],d=rightd,col='red',alpha=0.25)

	#plotcone()
	#points3d(cutcat[cutcat[,5]==1,1],cutcat[cutcat[,5]==1,2],cutcat[cutcat[,5]==1,3],size=5,col='red')

	return(cutcat)
}

#_________________________________________#
#_________________________________________#

# prepData
# Really quick function to take raw G3C catalogues and convert them into catalogues with 3D comoving Cartesian coordinates. Utilises CreateSample and limitCone. cat should either be a galaxy cat or a group cat. z is the redshift limit, c is the distance from the border. Medianz is a true/false flag that tells the code if it should take galaxies in groups and change their redshifts to the group's median redshift.

	prepData <- function(cat,field,z,c,appmaglim=19.4,mock=F, medianz=F){

	if(field=='G09'){rotangle=45;ralim=c(129,141);declim=c(-2,3);fieldn="09"}
	if(field=='G12'){rotangle=90;ralim=c(174,186);declim=c(-3,2);fieldn="12"}
	if(field=='G15'){rotangle=127.5;ralim=c(211.5,223.5);declim=c(-2,3);fieldn="15"}

	if(mock == T){groupcol = 8}else{groupcol = 9}

# If medianz is true, then we go through each galaxy in the catalogue. If it belongs in a group, it is assigned the redshift of the group. 
# WHAT IS THIS WHY DID I CODE THIS IT MAKES NO SENSE IGNORE

	if(medianz == T){
		for(i in 1:length(cat[,1])){
			if(cat[i,groupcol] > 0){cat[i,4] = G3CFoFGroup[which(G3CFoFGroup$GroupID==cat[i,groupcol]),'Zfof']}
		}
	cat(paste("Medianz applied","\n"))
	}

# Now calculate absolute magnitude limit for GAMA shallow fields, for given redshift.

	k = KEcorr(z)[1]
	DM = 5*log10(CosDist(z,H0=100,OmegaM=0.25,OmegaL=0.75)[,'LumDist']*1E5)
	M = appmaglim - DM - k

# First check to see if this is a groupcat or a galcat.

	if(colnames(cat)[1] == 'GroupID' | colnames(cat)[1] == 'HaloID'){type='group'}else{type='gal'}

# First chunk of code is for galaxies.

	if(type=='gal'){

# First do a magnitude cut by redshift. Calculating all the k-corrections takes ages so use approx.

		kfunction = approxfun(x=seq(0,0.6,by=0.0005),sapply(seq(0,0.6,by=0.0005),KEcorr)[1,])
		gal_k = sapply(cat[,4],kfunction)
		
		GalAbsMag = cat$Rpetro - cat$DM_100_25_75 - gal_k

		rows = which(cat$RA >= min(ralim) & cat$RA <= max(ralim) & cat[,4] <= z & GalAbsMag <= M)

		temp = rotate3d(sph2car(cat[rows,2],cat[rows,3],CosDist(cat[rows,4],H0=100,OmegaM=0.25,OmegaL=0.75)[,'CoDist']), (rotangle*pi/180),0,0,1)
		if(colnames(cat)[1]=='GalID'){temp = cbind(temp,GalAbsMag[rows],cat[rows,c(1,8)])}else{temp = cbind(temp,GalAbsMag[rows],cat[rows,c(1,9)])} 

# Check if mock or real galaxies. Depending on the result, catalogue will be GXX or MXX galxyz.

		if(colnames(cat)[1]=='GalID'){temp=cbind(temp,cat$Volume[rows]);colnames(temp)=c("X","Y","Z","AbsMag","GalID","HaloID","Volume")}else{colnames(temp)=c("X","Y","Z","AbsMag","CATAID","GroupID")}

	}

	if(type=='group'){

		rows = which(cat$IterCenRA >= min(ralim) & cat$IterCenRA <= max(ralim))
		tempcat = cat[rows,]

# Chuck the region into createSample to remove all galaxies and groups that cannot be seen past 0.212 (magnitude limited). Also, if this is a halo catalogue, the redshift column is differently labelled. It's always after the IterCenZ column though (halo or FoF mock) so just set that as the redshift column.

		

		if(mock==T){zcol = grep("IterCenZ",colnames(G3CMockFoFGroup))+1;zgroups = createSample(z,tempcat,G3CMockGal,appmaglim)}
		if(mock==F){zcol = grep("IterCenZ",colnames(G3CFoFGroup))+1;zgroups = createSample(z,tempcat,G3CGal,appmaglim)}
		keptgroups = tempcat[which(zgroups[,2]==0 & tempcat[,zcol] <= z),]
		cutgroups = limitCone(keptgroups,ralim,declim,c(0,z),c)

# Rotation is done in cutgroups, so ready to output!
		if(mock==T){cutgroups = cbind(cutgroups[order(cutgroups$GroupID),],keptgroups[order(keptgroups[,1]),'Volume'])}
		temp = cutgroups
		if(mock==T){colnames(temp) = c("X","Y","Z","GroupID","Cut","Volume")}

	}

	return(temp)

	}

#_________________________________________#
#_________________________________________#

# suppressMass
# This function works in conjunction with mass functions supplied by Steven Murray & Chris Power and is designed to artificially suppress the mass function of the GAMA simulation cones. It does this by calculating the predicted number density of a given group in different types of dark matter (using said mass functions) and uses the ratio of this predicted density to the current density to obtain a probability of that group existing in a different cosmology. It then generates a random number and discards the galaxy based on it.
# dmfunc should either be WDM1func or WDM01 func. Also need to make sure to have CDM func. Cat should be groups used to generate MST, groupcat should be whole catalogue (i.e. G3CMockHaloGroup).
# NOTE: In the HALO catalogues for making trees, the Halo ID is called GroupID!

suppressMass <- function(dmfunc,cat,groupcat){

# Select out groups.

	groups = groupcat[which(groupcat$HaloID %in% cat$GroupID),]
	keptgroups = {}

# Loop over each group and work out density in CDM and target DM model.

	for(i in 1:dim(groups)[1]){

		CDM = CDMfunc(log10(groups[i,'HaloMass']))
		DM = dmfunc(log10(groups[i,'HaloMass']))
		ratio = CDM/DM

# ratio is the value that the random needs to be over. If it's greater than 1, keep the halo by default.

		if(ratio <= 1){

			if(runif(1) <= ratio){keptgroups = rbind(keptgroups,groups[i,])}else{next}

		}else{keptgroups = rbind(keptgroups,groups[i,])}
	}

	outcat = cat[which(cat$GroupID %in% keptgroups$HaloID),]
	return(outcat)
}

#_________________________________________#
#_________________________________________#

## nodePrune ##
# Runs a k-node prune and edge cut on a minimal spanning tree generated using the mst() function in the nnclust package. Spits out a reduced MST as a list with from, to and dist. cat should be a data frame containing X, Y, Z coordinates of particles used to make the tree.

nodePrune <- function(tree,k,b,cat,plot=F){

# Begin by calculating the number of links each node has. Do this by looping through tree$from and counting the number of times each element appears.

	links = {};for(i in 1:length(tree$from)){links=rbind(links,length(tree$from[which(tree$from==tree$from[i])]))}
	
# Can now reduce MST by using which(). Remember, mst() outputs squared distances so must square distance limit.

	reduced_to = tree$to[which(links <= k & tree$dist <= (b^2))]
	reduced_from = tree$from[which(links <= k & tree$dist <= (b^2))]
	reduced_dist = tree$dist[which(links <= k & tree$dist <= (b^2))]
	
	reduced_mst = list(from = reduced_from, to = reduced_to, dist = reduced_dist)
	
# Optional; plotting.

	if(plot){
		aplot(cat$X[unique(c(reduced_mst$from,reduced_mst$to))],cat$Y[unique(c(reduced_mst$from,reduced_mst$to))],asp=1,pch='.',col='blue',xlab=expression(paste("x (h"^{"-1"}, "Mpc)",sep="")),ylab=expression(paste("y (h"^{"-1"}, "Mpc)",sep="")),cex=0.5)
		
	segments(cat$X[reduced_mst$from],cat$Y[reduced_mst$from],cat$X[reduced_mst$to],cat$Y[reduced_mst$to],col=hsv(h=1,v=0,alpha=0.5))
	
	}

	return(reduced_mst)

}

#_________________________________________#
#_________________________________________#

## makePrimCat ##
# This function generates the primary filament catalogue, consisting of groups linked together using mst, nodePrune and groupit. Needs the MST, k and b, and the group catalogue the tree is made from. Catalogue must contain a GroupID column. Gals must contain CATAIDs as well as 3D cartesian coords. Walkargs is a three-element vector containing logical arguments (TRUE or FALSE) for bydist, dodist and a 3D cartesian data frame for xyz; all to be fed into makenetwork if necessary.

# b is the maximum edge length in Mpc, r is the scooping distance for picking up galaxies. Prefix is just a tag that goes in front of the treeID to distinguish between different cones/mock volumes; must be a number.

makePrimCat <- function(tree,k=Inf,b,cat,gals,r,q,filByDist = F,vol=0,halo=F,prefix=""){

	filaments = list(links={},network={})
	groupcat = {}
	filcat = {}
	treecat = {}
	galcat = {}

# If running on mocks, reduce galaxy data frame down to just the volume being considered.

	if(vol > 0){gals = gals[which(gals$Volume==vol),]}

# First run nodePrune and groupit. This will generate a pruned tree and group objects together.

	tempprune = nodePrune(tree,k=k,b=b,cat=cat,plot=F)
	tempgroup = groupit(cbind(c(tempprune$to,tempprune$from),c(tempprune$from,tempprune$to)))

# Generate a list where each entry is an individual filament.

	for(i in 1:length(unique(tempgroup[,2]))){

# Generate filament catalogue containing groupID and filament ID. Catalogue is a list where each element contains from and to links for a filament. To get just the filament, do a unique search on a vector where both columns have been appended together.

		links = cbind(tempprune$from[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]],tempprune$to[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]],tempprune$dist[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]])
		links = cbind(links[,1],links[,2],sqrt(links[,3]))
		
		# Run makenetwork for each filament and append its results, if there are 2 or more groups. If there are 2 groups, manually construct network.

		if(length(links[,1])>1){network = makenetwork(links=links[,c(1,2)],bydist=filByDist,dodist=T,xyz=cbind(unique(c(links[,1], links[, 2])),cat[unique(c(links[,1], links[, 2])), c("X", "Y", "Z")]))}
		else{
			fb = list(fil=c(links[1,1],links[1,2]),links=links)
			backbone = c(links[1,1],links[1,2])
			Nbb = 2
			bblinks = links
			bbdist = rdist(cat[links[1,1],c("X","Y","Z")],cat[links[1,2],c("X","Y","Z")])
			bblist = bblist = list(backbone = backbone, Nbb = Nbb, bblinks = bblinks, bbdist = bbdist)
			branchlist = c(1,1)
			network=list(fils=fb,bblist=bblist,network=branchlist)}

		filaments$links[[i]]=links;filaments$network[[i]] = network
	
	}
	
# Now time to go through each sub-tree and tag them up if there are groups that are in the cut region.

	allcutgroups = {}

	for(i in 1:length(unique(tempgroup[,2]))){
		
# Begin by identifying which groups in the filament are cut and find which ones to keep.

		cutgroups = which(cat$GroupID %in% cat[unique(c(filaments$links[[i]][,1],filaments$links[[i]][,2])),'GroupID'])[which(cat[which(cat$GroupID %in% cat[unique(c(filaments$links[[i]][,1],filaments$links[[i]][,2])),'GroupID']),'Cut']==1)]

		if(length(cutgroups)>0){allcutgroups = append(allcutgroups,as.vector(cutgroups))}

		#keeplinks = filaments$links[[i]][!filaments$links[[i]][,1] %in% cutgroups & !filaments$links[[i]][,2] %in% cutgroups,c(1,2,3)]

	}

# Now run through the whole data set and make a number of catalogues. Also add columns called Cut to show if in cut region or not. 1 for cut, 0 for included.
# TreeCat contains a treeID, the number of branches in each tree, the number of groups, and if it has any cut groups.
# FilCat contains individual IDs for each branch in the network of each tree, the tree it's in, the number of groups, and if any are cut.
# GroupCat contains a GroupID, a branchID, the filament ID, and the tree ID, and if it is cut.

	print("Generating filaments...")

	for(i in 1:length(filaments$links)){
		tempdist = {}

# Clunky bit of code to format treeID. Basically add 0s to the front of each treeID such that they have as many digits as the total number of trees.

		treeID = as.numeric(paste(prefix, formatC(i, flag = 0, width = 4),sep=""))
		Ntgroup = length(unique(c(filaments$links[[i]][,1],filaments$links[[i]][,2])))

# Check here if tree has more than 2 groups. If not, it can only have 1 branch.

		if(Ntgroup==2){Nbranch=1}else{Nbranch = length(filaments$network[[i]]$fils)}
		Ntcut = length(which(unique(c(filaments$links[[i]][, 1], filaments$links[[i]][,2])) %in% allcutgroups==T))	
		
# Also gather up the total length of all the links in the tree. Do this by quickly summing over each branch.
		if(Ntgroup >2){
			for(k in 1:length(filaments$network[[i]]$fils)){

				filDist = sum(rbind(filaments$network[[i]]$fils[[k]]$links[,'linkdist']))
				tempdist = c(tempdist, filDist)
			}
		}else{tempdist = filaments$network[[i]]$fils$links[,3]}

		tlength = sum(rbind(tempdist))

# Grab the backbone length too.

		bblength = filaments$network[[i]]$bblist$bbdist
		treecat = rbind(treecat,c(treeID, Nbranch, Ntgroup, tlength, bblength, Ntcut))

# Now construct filcat.
	
		for(j in 1:Nbranch){

			temp = formatC(j, flag=0, width=4)
			filID = as.numeric(paste(treeID,temp,sep=""))

# Again, check for length of tree.

			if(Ntgroup==2){Nfgroup = 2;order=1;Nfcut = Ntcut;flength = filaments$network[[i]]$fils$links[,3]}else{
			Nfgroup = filaments$network[[i]]$fils[[j]]$N
			order = filaments$network[[i]]$network[which(filaments$network[[i]]$network[,1]==j),2]
			Nfcut = length(which(filaments$network[[i]]$fils[[j]]$fil %in% allcutgroups==T))
			flength = sum(rbind(filaments$network[[i]]$fils[[j]]$links[,'linkdist']))}

			filcat = rbind(filcat,c(filID,treeID,order, Nfgroup, flength, Nfcut))

		}

# Now construct groupcat. For the current tree, list all groups present in it. Then skim through each branch to identify the branch it's in to reconstruct the filID, and also pull out the groupID. Also check to see if it's cut.
	
		treegroups = unique(c(filaments$links[[i]][,1],filaments$links[[i]][,2]))
		for(j in 1:length(treegroups)){
			if(Ntgroup==2){gfilID = as.numeric(paste(treeID,formatC(1, flag=0, width=4),sep=''));gorder = 1;groupID = cat$GroupID[treegroups[j]]}else{
			k = 1
			groupID = cat$GroupID[treegroups[j]]
		
			while(any(filaments$network[[i]]$fils[[k]]$fil==treegroups[j])==F){k=k+1}
			
			gfilID = as.numeric(paste(treeID,formatC(k, flag=0, width=4),sep=""))
			gorder = as.numeric(filaments$network[[i]]$network[which(filaments$network[[i]]$network[,1]==k),2])}

			gcheckifcut = as.numeric(cat$Cut[treegroups[j]])
			groupinfo = as.numeric(G3CFoFGroup[which(cat$GroupID[treegroups[j]] == G3CFoFGroup$GroupID),c("IterCenRA","IterCenDec","Zfof")])
			groupcat = rbind(groupcat,c(groupID,groupinfo,as.numeric(cat[treegroups[j],c("X","Y","Z")]),gfilID,gorder,treeID,gcheckifcut))

		}

	}

# filaments$cat is the final list, containing GroupID-FilID.

		colnames(treecat) = c("TreeID","Nbranch","Ngroup", "Length", "BBlength", "Cut")
		colnames(filcat) = c("FilID","TreeID","Order","Ngroup", "Length", "Cut")
		colnames(groupcat) = c("GroupID", "RA", "Dec", "Zfof", "X", "Y", "Z","FilID","Order","TreeID","Cut")
		
		filaments$treecat = treecat
		filaments$filcat = filcat
		filaments$groupcat = groupcat
		filaments$cutgroups = allcutgroups

# Output the MST too, but make sure it contains GroupIDs and not just numbers.

		temptree = cbind(as.integer(cat$GroupID[tempprune$from]),as.integer(cat$GroupID[tempprune$to]),tempprune$dist)
		filaments$mst = cbind(paste(prefix,seq(1:dim(temptree)[1]),sep=""),temptree)
		colnames(filaments$mst) = c("LinkID","From","To","Dist")

# Next, construct galcat. For each tree, run scooper and put together the results. Only do this if a galcat is provided. If we are using mocks, then vol > 0 and we need to make sure we only select galaxies from the relevant volume.

# Now - any galaxy that belongs to any of these groups is immediately associated with that filament/branch. Their distance is the projected distance to the group centre.

print(paste("Generating galaxies..."))

	galsingroups = {}
	discardrows = {}

	for(i in 1:length(gals[,1])){

		if(vol==0){
			groupid = G3CGal[which(G3CGal$CATAID == gals$CATAID[i]),'GroupID']
		}else{
			if(halo==F){groupid = G3CMockGal[which(G3CMockGal$GalID == gals[i,'GalID']),'GroupID']}
			if(halo==T){groupid = G3CMockGal[which(G3CMockGal$GalID == gals[i,'GalID']),'HaloID']}
		}

		if(groupid == 0 | length(which(filaments$groupcat[,'GroupID']==groupid)) == 0){
			next}
		else{

			groupcen = filaments$groupcat[which(filaments$groupcat[,'GroupID']==groupid),c("X","Y","Z")]
			galcen = gals[i,c("X","Y","Z")]

			# Get projected distance.

			projdist = as.numeric(sqrt((galcen[1]-groupcen[1])^2 + (galcen[2]-groupcen[2])^2 + (galcen[3]-groupcen[3])^2))

			# Add this galaxy to list, and the row from the rows to be discarded.
			if(projdist <= r){
				if(vol == 0){
					galsingroups = rbind(galsingroups, c(gals[i,'CATAID'],as.numeric(galcen),as.numeric(projdist),filaments$groupcat[which(filaments$groupcat[,'GroupID']==groupid),c("FilID","TreeID")]))
					}else{
					galsingroups = rbind(galsingroups, c(gals[i,'GalID'],as.numeric(galcen),as.numeric(projdist),filaments$groupcat[which(filaments$groupcat[,'GroupID']==groupid),c("FilID","TreeID")]))
					}

				discardrows = c(discardrows,i)
			}
		}
	}

	if(vol==0){colnames(galsingroups) = c("CATAID","X","Y","Z","d","FilID","TreeID")}
	if(vol >0){colnames(galsingroups) = c("GalID","X","Y","Z","d","FilID","TreeID")}

	# Remove rows from gals.

	gals = gals[-discardrows,]

	if(vol==0){
		if(length(gals) > 1){

			for(i in 1:length(filaments$links)){
				
				temp = scooper(filaments$treecat[i,'TreeID'], r, filaments, cat, gals)
				if(length(temp)>0){galcat = rbind(galcat,cbind(temp,as.numeric(filaments$treecat[i,'TreeID'])))}
			}
						
			colnames(galcat) = c("CATAID","X","Y","Z","d","FilID","TreeID")
			galcat = rbind(galcat,galsingroups)
			
		}
	}

	if(vol > 0){
		if(length(gals[which(gals$Volume==vol),]) > 1){

			for(i in 1:length(filaments$links)){
				temp = scooper(filaments$treecat[i,'TreeID'], r, filaments, cat, gals[which(gals$Volume==vol),])
				if(length(temp)>0){galcat = rbind(galcat,cbind(temp,as.numeric(filaments$treecat[i,'TreeID'])))}
			}
			
			colnames(galcat) = c("GalID","X","Y","Z","d","FilID","TreeID")
			galcat = rbind(galcat,galsingroups)
			
		}
	}

# Go through galcat and check for duplicates. Not sure yet why this is happening so this is a temporary fix. In cases of duplicates, select option with the lesser value of d.

	duplicates = {}
	temp = cbind(as.numeric(names(table(galcat[,1]))),as.numeric(table(galcat[,1])))
	temp = temp[which(temp[,2]>1),]

	# Only proceed if there actually ARE duplicates. Unlikely that there wouldn't be, but it happens.

	if(dim(temp)[1]>0){

		for(i in 1:dim(temp)[1]){
			duplicates = galcat[which(galcat[,1] == temp[i,1]),]
			
			#Identify row(s) to be removed.
			killrow = duplicates[-which.min(duplicates$d),]

			#Find row in galcat and OBLITERATE!!
			if(dim(killrow)[1]==1){
				galcat = galcat[-which(galcat$d == killrow[,'d'] & galcat[,1] == killrow[,1]),]
			}
			if(dim(killrow)[1]>1){
				for(j in 1:dim(killrow)[1]){
					galcat = galcat[-which(galcat$d == killrow[j,'d'] & galcat[,1] == killrow[j,1]),]
				}	
			}
		}
	}
	
	filaments$galcat = galcat

# Now construct the tendril catalogue. Gather up unscooped galaxies and make an MST out of them.

	print("Generating tendrils...")
	
	if(vol == 0){unscoopedgals = gals[which(!gals[,'CATAID'] %in% galcat[,1]),]}
	if(vol > 0){unscoopedgals = gals[which(!gals[which(gals$Volume==vol),'GalID'] %in% galcat[,'GalID']),]}
	tree = mst(as.matrix(unscoopedgals[,c("X","Y","Z")]))

# Now prune tree and regroup. q determines tendril MST cut off.

	tempprune = nodePrune(tree,k=k,b=q,cat=unscoopedgals,plot=F)
	tempgroup = groupit(cbind(c(tempprune$to,tempprune$from),c(tempprune$from,tempprune$to)))

# Now generate a tendril catalogue. For each tendril specify number of galaxies, total length, and give it a unique ID.

	tendrilcat = {}
	tendrilgals = {}
	voidgals = {}

	for(i in 1:length(unique(tempgroup[,2]))){

		links = cbind(tempprune$from[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]],tempprune$to[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]],tempprune$dist[tempprune$from %in% tempgroup[which(tempgroup[,'filID']==i),1]])
		links = cbind(links[,1],links[,2],sqrt(links[,3]))

		tendrilLength = sum(links[,3])
		tendrilN = length(unique(c(links[,1],links[,2])))
		tendrilID = as.numeric(paste(prefix, formatC(i, flag = 0, width = 4),sep=""))

		tendrilcat = rbind(tendrilcat, c(tendrilID,tendrilN,tendrilLength))
		tendrilgals = rbind(tendrilgals, cbind(unscoopedgals[unique(c(links[,1],links[,2])),], tendrilID))
				
	}

	colnames(tendrilcat) = c("TendrilID","Ngal","Length")
	if(vol==0){tendrilgals = tendrilgals[,c(1,2,3,4,5,7)];colnames(tendrilgals) = c("X","Y","Z","AbsMag","CATAID","TendrilID")}
	if(vol>0){tendrilgals = tendrilgals[,c(1,2,3,4,5,8)];colnames(tendrilgals) = c("X","Y","Z","AbsMag","GalID","TendrilID")}

	# Any group that is not in the MST generated to make the tendrils is classified as a 'void' galaxy.
	
	voidgals = unscoopedgals[-unique(tempgroup[,1]),]

	# Remove all void galaxies that are in groups.

	voidgals = voidgals[which(voidgals$GroupID == 0),]

	filaments$tendrilcat = tendrilcat
	filaments$tendrilgals = tendrilgals

	# Output tendril MST as well, this time with CATAIDs.
	if(vol==0){
		tendriltemptree = cbind(as.integer(unscoopedgals$CATAID[tempprune$from]),as.integer(unscoopedgals$CATAID[tempprune$to]),tempprune$	dist)
		filaments$tmst = cbind(paste(prefix,seq(1:dim(tendriltemptree)[1]),sep=""),tendriltemptree)
		
	}
	if(vol > 0){
		tendriltemptree = cbind(as.integer(unscoopedgals$GalID[tempprune$from]),as.integer(unscoopedgals$GalID[tempprune$to]),tempprune$dist)
		filaments$tmst = cbind(paste(prefix,seq(1:dim(tendriltemptree)[1]),sep=""),tendriltemptree)
	}

	colnames(filaments$tmst) = c("LinkID","From","To","Dist")

	filaments$voidgals = voidgals

	return(filaments)
}

#_________________________________________#
#_________________________________________#

## scooper ##
# This function takes a filament and scoops up all galaxies around each member of the filament, out to a distance b. By searching for galaxies at a distance r from a line connecting each node, it scoops up galaxies along the line and at spheres around each point. Keeps track of orthogonal distance from galaxy to nearest node or edge. Read in fil from the filaments that is spat out by makePrimCat.

# Cat should be the output of makePrimCat. groupxyz should be the xyz catalogue used to generate the tree for that region.

scooper <- function(treeID, r, cat, groupxyz, galcat){

	# Start by identifying which groups are in the filament.

	treegroups = cat$groupcat[which(cat$groupcat[,'TreeID']==treeID),]
	treegroups = as.data.frame(treegroups)

	# Now go through links, defining end points as the group at each end of the link. This is all that is needed.
	# Usual count-in-spheres around the end points, and distance between a line and a point formula for the cylinder.

	cylgals = {}
	for(i in 1:dim(cat$links[[which(cat$treecat[,'TreeID']==treeID)]])[1]){

		sphcen1 = treegroups[treegroups$GroupID == groupxyz[cat$links[[which(cat$treecat[,'TreeID']==treeID)]][i,1],'GroupID'],]
		sphcen2 = treegroups[treegroups$GroupID == groupxyz[cat$links[[which(cat$treecat[,'TreeID']==treeID)]][i,2],'GroupID'],]

		# Define a 'box' around the link, and use this as a subset when searching for galaxies. All galaxies beyond these limits are excluded.

		xmax = max(as.numeric(sphcen1$X,sphcen2$X))+r; xmin = min(as.numeric(sphcen1$X,sphcen2$X))-r
		ymax = max(as.numeric(sphcen1$Y,sphcen2$Y))+r; ymin = min(as.numeric(sphcen1$Y,sphcen2$Y))-r
		zmax = max(as.numeric(sphcen1$Z,sphcen2$Z))+r; zmin = min(as.numeric(sphcen1$Z,sphcen2$Z))-r

		galsel = which(galcat$X >= xmin & galcat$X <= xmax & galcat$Y >= ymin & galcat$Y <= ymax & galcat$Z >= zmin & galcat$Z <= zmax)

		# If galsel is empty, that means that all galaxies near this filament have already been associated with a filament. Exit.
		
		if(length(galsel)==0){next}

		# Now select galaxies within a distance r of the line connecting both spheres.

		# Start by defining the vector the two points are on (P1 and P2). For this need one point, P1, and the vector joining the two, u. The line is called PP'

		P1 = matrix(c(sphcen1$X,sphcen1$Y,sphcen1$Z),ncol=3)
		P2 = matrix(c(sphcen2$X,sphcen2$Y,sphcen2$Z),ncol=3)
		u = P1 - P2

		# Now loop through galsel and work out distance for each point. Each galaxy is Q.
		
		d = {}
		for(i in 1:length(galsel)){
			
			Q = matrix(c(galcat$X[galsel[i]],galcat$Y[galsel[i]],galcat$Z[galsel[i]]),ncol=3)
			v = P1 - Q
			
			c1 = (as.vector(u)%*%as.vector(v))/sqrt(sum(u^2))^2

			# Want the projection of Q on the line, which is P

			Q_L = P1 + c1*as.vector(u)

			dQQ_L = sqrt((Q[1]-Q_L[1])^2+(Q[2]-Q_L[2])^2+(Q[3]-Q_L[3])^2)
			dQP1 = sqrt((Q[1]-P1[1])^2+(Q[2]-P1[2])^2+(Q[3]-P1[3])^2)
			dQP2 = sqrt((Q[1]-P2[1])^2+(Q[2]-P2[2])^2+(Q[3]-P2[3])^2)

			d = c(d,min(c(dQQ_L,dQP1,dQP2)))
			}
		
		# normalLength = sqrt((sphcen2$X[[1]]-sphcen1$X[[1]])^2 + (sphcen2$Y[[1]] - sphcen1$Y[[1]])^2)
		# d = abs((galcat$X[galsel]-sphcen1$X[[1]])*(sphcen2$Y[[1]]-sphcen1$Y[[1]]) - (galcat$Y[galsel]-sphcen1$Y[[1]])*(sphcen2$X[[1]]-sphcen1$X[[1]]))/normalLength

		IDcolname = colnames(galcat)[which(colnames(galcat) == 'CATAID' | colnames(galcat) == 'GalID')]

		if(length(galsel)!=0 & length(which(d<=r)) > 0){
			cylgals = rbind(cylgals,cbind(galcat[galsel[which(d <= r)],c(IDcolname,"X","Y","Z")],as.numeric(d[which(d <= r)]),as.numeric(sphcen2$FilID)))
		}

	}
	
	# Only do next bit if any galaxies are actually scooped up.

	while(length(cylgals)>0){
		colnames(cylgals) = c("CATAID","X","Y","Z","d","FilID")

		# All galaxies now scooped up. Next step is to go through them all and detect degeneracies (galaxies scooped up by more than one branch. In these cases, take the nearest cylinder.)

		temp = table(cylgals$CATAID)
		degens = cbind(as.numeric(names(temp)),as.numeric(temp))
		nodegens = degens[which(degens[,2]==1),]
		degens = degens[which(degens[,2] > 1),]
		newcylgals = {}

		if(length(degens)>2){

			for(i in 1:length(degens[,1])){

				tempgals = cylgals[which(cylgals$CATAID == degens[i,1]),]
				keepgal = tempgals[which.min(tempgals$d),]

				# Stick galaxies to keep into a new data frame.

				newcylgals = rbind(newcylgals,keepgal)

			}

		}else{

			tempgals = cylgals[which(cylgals$CATAID == degens[1]),]
			keepgal = tempgals[which.min(tempgals$d),]
			newcylgals = rbind(newcylgals,keepgal)
		}

		# Now combine newcylgals with all cylgal galaxies that have no degeneracies.

		if(length(nodegens)>2){newcylgals = rbind(newcylgals,cylgals[cylgals$CATAID %in% nodegens[,1],])}else{
			newcylgals = rbind(newcylgals,cylgals[cylgals$CATAID %in% nodegens[1],])}
		
		newcylgals = newcylgals[order(newcylgals$FilID),]

		# All done. Output!
		return(newcylgals)
	}
}

#_________________________________________#
#_________________________________________#

## plotFil ##
# Quick and dirty function to plot a primary filament made with makePrimCat. Filament galaxies are shown in blue, unassociated groups in red, and galaxies in black (all galaxies shown, not just ungrouped ones). Edges are drawn to show the backbone as thick and red, with other branches shown in decreasing thickness and more yellow.

plotFil <- function(filID, filcat, groupcat, galcat, dim=2, text=F, bydist=F, dodist=F, xyz={}){

	filIDOrig = filID
	filID = which(filcat$treecat[,'TreeID']==filID)
	links = filcat$links[[filID]][,c(1,2)]
	branches = filcat$network[[filID]]
	#cols = hsv(h=1,s=seq(1,0.05,by=-1/max(branches$network[,2])))
	#cols = rainbow(n=max(branches$network[,2])+2,alpha=1)
	cols = rev(brewer.pal(5,"Paired"))
	
	central_tree = filcat$network[[filID]]$walk[which.max(filcat$network[[filID]]$walk[,2]),1]

	if(dim==2){
		aplot(groupcat$X[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Y[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],pch='.',col='blue',asp=1,cex=3,xlab=expression(paste("RA [h"^{"-1"}, " Mpc]")),ylab=expression(paste("z [h"^{"-1"}, " Mpc]")))

		points(groupcat$X[-unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Y[-unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],pch='.',col='green4',asp=1,cex=3)

		points(galcat$X,galcat$Y,pch='.',col=hsv(h=1,v=0,alpha=0.5),asp=1)

		for(i in length(branches$network[,2]):3){

			templinks = links[which(links[,1] %in% branches$fils[[branches$network[i,1]]]$fil | links[,2] %in% branches$fils[[branches$network[i,1]]]$fil),]
			if(length(templinks)==2){templinks=matrix(data=templinks,nrow=1,ncol=2)}
			segments(groupcat$X[templinks[,1]],groupcat$Y[templinks[,1]],groupcat$X[templinks[,2]],groupcat$Y[templinks[,2]],lwd=3,col=cols[branches$network[i,2]])
		}

		segments(groupcat$X[branches$bb$bblinks[,1]],groupcat$Y[branches$bb$bblinks[,1]],groupcat$X[branches$bb$bblinks[,2]],groupcat$Y[branches$bb$bblinks[,2]],lwd=5,col=cols[1])

		points(groupcat$X[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Y[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],pch='.',col='blue',asp=1,cex=3)

		label(pos='topleft',lab=paste("N = ",dim(filcat$links[[filID]])[1]+1,sep=""),bty='b')
	}

	if(dim==3){
		require(rgl)

		plot3d(groupcat$X[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Y[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Z[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],size=5.5,col='blue',xlab='',ylab='',zlab='',axes=F)
		

		if(text){text3d(x=groupcat$X[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],y=groupcat$Y[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],z=groupcat$Z[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],texts=unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2])),adj=c(0,0),cex=0.75)}

		if(length(unique(rbind(branches$network)[,2])) > 1){

			for(i in length(branches$network[,2]):3){

				templinks = links[which(links[,1] %in% branches$fils[[branches$network[i,1]]]$fil | links[,2] %in% branches$fils[[branches$network[i,1]]]$fil),]

				if(length(templinks)==2){templinks=matrix(data=templinks,nrow=1,ncol=2)}

				for(j in 1:length(templinks[,1])){

					#lines3d(x=c(groupcat$X[templinks[j,1]],groupcat$X[templinks[j,2]]),y=c(groupcat$Y[templinks[j,1]],groupcat$Y[templinks[j,2]]),z=c(groupcat$Z[templinks[j,1]],groupcat$Z[templinks[j,2]]),lwd=150,col=cols[branches$network[i,2]])

					mat = matrix(data=cbind(c(groupcat$X[templinks[j,1]],groupcat$X[templinks[j,2]]),c(groupcat$Y[templinks[j,1]],groupcat$Y[templinks[j,2]]),c(groupcat$Z[templinks[j,1]],groupcat$Z[templinks[j,2]])),nrow=2,ncol=3,byrow=F)

					temp = cylinder3d(center=mat,radius=0.2,sides=20)
					shade3d(temp,col=cols[branches$network[i,2]])
				}

			}
		}	
		
		#temp = {}
		for(i in 1:length(branches$bb$bblinks[,1])){
			#lines3d(x=c(groupcat$X[branches$bb$bblinks[i,1]],groupcat$X[branches$bb$bblinks[i,2]]),y=c(groupcat$Y[branches$bb$bblinks[i,1]],groupcat$Y[branches$bb$bblinks[i,2]]),z=c(groupcat$Z[branches$bb$bblinks[i,1]],groupcat$Z[branches$bb$bblinks[i,2]]),lwd=100,col=cols[1])

			mat = matrix(data=cbind(c(groupcat$X[branches$bb$bblinks[i,1]],groupcat$X[branches$bb$bblinks[i,2]]),c(groupcat$Y[branches$bb$bblinks[i,1]],groupcat$Y[branches$bb$bblinks[i,2]]),c(groupcat$Z[branches$bb$bblinks[i,1]],groupcat$Z[branches$bb$bblinks[i,2]])),nrow=2,ncol=3,byrow=F)

			temp = cylinder3d(center=mat,radius=0.3,sides=20)
			shade3d(temp,col=cols[1])
			#temp = rbind(temp,c(groupcat$X[branches$bb$bblinks[i,1]],groupcat$X[branches$bb$bblinks[i,2]],groupcat$Y[branches$bb$bblinks[i,1]],groupcat$Y[branches$bb$bblinks[i,2]],groupcat$Z[branches$bb$bblinks[i,1]],groupcat$Z[branches$bb$bblinks[i,2]]))
		}
		
		#colnames(temp)=c("X1","X2","Y1","Y2","Z1","Z2")
		#write.csv(temp,paste("~/Desktop/bblinkcats/",filIDOrig,"_BBlinks.csv",sep=""),row.names=F)

		spheres3d(groupcat$X[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Y[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],groupcat$Z[unique(append(filcat$links[[filID]][,1],filcat$links[[filID]][,2]))],radius=0.5,col='blue')
		spheres3d(groupcat$X[central_tree],groupcat$Y[central_tree],groupcat$Z[central_tree],radius=1,col=metCols[6])

		#cutgroups = groupcat[unique(append(filcat$links[[filID]][, 1], filcat$links[[filID]][,2]))[which(groupcat[unique(append(filcat$links[[filID]][, 1], filcat$links[[filID]][,2])),'Cut']==1)],]
		#if( dim(cutgroups)[1]>0 ){spheres3d(cutgroups$X,cutgroups$Y,cutgroups$Z,col='black',radius=0.35)}

	}
}

#_________________________________________#
#_________________________________________#

plotSliceZoom  <- function(ralim,declim,zlim,region,vol=0,labels=F,mst=T,pgals=T,tendrils=T,pvoids=T,lowZ=F){

	#make sure ralim,declim, and zlim are all min-max ranges

	#in case you want to do mocks, use this code

	if(vol>0){
		if(region=="G09"){cat = get(paste("M09b5.75c1v",vol,"cat",sep=""));groups=M09c1xyz[which(M09c1xyz$Volume==vol),];gals=M09zgal[which(M09zgal$Volume==vol),];rotangle=45}
		if(region=="G12"){cat = get(paste("M12b5.75c1v",vol,"cat",sep=""));groups=M12c1xyz[which(M12c1xyz$Volume==vol),];gals=M12zgal[which(M12zgal$Volume==vol),];rotangle=90}
		if(region=="G15"){cat = get(paste("M15b5.75c1v",vol,"cat",sep=""));groups=M15c1xyz[which(M15c1xyz$Volume==vol),];gals=M15zgal[which(M15zgal$Volume==vol),];rotangle=127.5}

		yrange = c(CosDist(min(zlim),100,0.25,0.75)[,'CoDist'],CosDist(max(zlim),100,0.25,0.75)[,'CoDist'])

		filgroupIDs = cat$groupcat[which(cat$groupcat[,'GroupID'] %in% G3CMockFoFGroup$GroupID[which(G3CMockFoFGroup$IterCenRA >= min(ralim) & G3CMockFoFGroup$IterCenRA <= max(ralim) & G3CMockFoFGroup$IterCenDec >= min(declim) & G3CMockFoFGroup$IterCenDec <= max(declim) & G3CMockFoFGroup$Zfof >= min(zlim) & G3CMockFoFGroup$Zfof <= max(zlim))]),'GroupID'] 
		filgalIDs = cat$galcat[which(cat$galcat[,'GalID'] %in% G3CMockGal$GalID[which(G3CMockGal$RA >= min(ralim) & G3CMockGal$RA <= max(ralim) & G3CMockGal$DEC >= min(declim) & G3CMockGal$DEC <= max(declim) & G3CMockGal$Zspec >= min(zlim) & G3CMockGal$Zspec <= max(zlim))]),'GalID'] 
		tendrilIDs = cat$tendrilgals[which(cat$tendrilgals[,'GalID'] %in% G3CMockGal$GalID[which(G3CMockGal$RA >= min(ralim) & G3CMockGal$RA <= max(ralim) & G3CMockGal$DEC >= min(declim) & G3CMockGal$DEC <= max(declim) & G3CMockGal$Zspec >= min(zlim) & G3CMockGal$Zspec <= max(zlim))]),'GalID'] 
		voidIDs = cat$voidgals[which(cat$voidgals[,'GalID'] %in% G3CMockGal$GalID[which(G3CMockGal$RA >= min(ralim) & G3CMockGal$RA <= max(ralim) & G3CMockGal$DEC >= min(declim) & G3CMockGal$DEC <= max(declim) & G3CMockGal$Zspec >= min(zlim) & G3CMockGal$Zspec <= max(zlim))]),'GalID'] 

		#plot all groups in filaments

		plot(-cat$groupcat[which(cat$groupcat[,'GroupID'] %in% filgroupIDs),'X'],cat$groupcat[which(cat$groupcat[,'GroupID'] %in% filgroupIDs),'Y'],cex=0.5,asp=1,axes=F,xlab="",ylab="",ylim=c(min(yrange)-10,max(yrange)+20),col='white',pch='.')
		if(pgals==T){points(-cat$galcat[which(cat$galcat[,'GalID'] %in% filgalIDs),'X'],cat$galcat[which(cat$galcat[,'GalID'] %in% filgalIDs),'Y'],pch='.',col='blue',cex=4)}
		if(tendrils==T){points(-cat$tendrilgals[which(cat$tendrilgals[,'GalID'] %in% tendrilIDs),'X'],cat$tendrilgals[which(cat$tendrilgals[,'GalID'] %in% tendrilIDs),'Y'],pch='.',col='green4',cex=4)}
		if(pvoids==T){points(-cat$voidgals[which(cat$voidgals[,'GalID'] %in% voidIDs),'X'],cat$voidgals[which(cat$voidgals[,'GalID'] %in% voidIDs),'Y'],pch='.',col='red',cex=4)}
		overlaygrid(ralim=ralim,declim=declim,zlim=zlim,text=c(1,2,3),textsize=0.75,radeclabels=labels)
	}else{

		#for real data
		if(lowZ == F){
			if(region=="G09"){cat = G09b5.75c1cat;groups=G09c1xyz;gals=G09zgal;rotangle=45}
			if(region=="G12"){cat = G12b5.75c1cat;groups=G12c1xyz;gals=G12zgal;rotangle=90}
			if(region=="G15"){cat = G15b5.75c1cat;groups=G15c1xyz;gals=G15zgal;rotangle=127.5}
		}
		if(lowZ == T){
			if(region=="G09"){cat = G09b3.7cat;groups=G09xyzLowZ;gals=G09galLowZ;rotangle=45}
			if(region=="G12"){cat = G12b3.7cat;groups=G12xyzLowZ;gals=G12galLowZ;rotangle=90}
			if(region=="G15"){cat = G15b3.7cat;groups=G15xyzLowZ;gals=G15galLowZ;rotangle=127.5}
		}

		yrange = c(CosDist(min(zlim),100,0.25,0.75)[,'CoDist'],CosDist(max(zlim),100,0.25,0.75)[,'CoDist'])

		filgroupIDs = cat$groupcat[which(cat$groupcat[,'GroupID'] %in% G3CFoFGroup$GroupID[which(G3CFoFGroup$IterCenRA >= min(ralim) & G3CFoFGroup$IterCenRA <= max(ralim) & G3CFoFGroup$IterCenDec >= min(declim) & G3CFoFGroup$IterCenDec <= max(declim) & G3CFoFGroup$Zfof >= min(zlim) & G3CFoFGroup$Zfof <= max(zlim))]),'GroupID']
		filgalIDs =  cat$galcat[which(cat$galcat[,'CATAID'] %in% G3CGal$CATAID[which(G3CGal$RA >= min(ralim) & G3CGal$RA <= max(ralim) & G3CGal$Dec >= min(declim) & G3CGal$Dec <= max(declim) & G3CGal$Z >= min(zlim) & G3CGal$Z <= max(zlim))]),'CATAID'] 
		tendrilIDs = cat$tendrilgals[which(cat$tendrilgals[,'CATAID'] %in% G3CGal$CATAID[which(G3CGal$RA >= min(ralim) & G3CGal$RA <= max(ralim) & G3CGal$Dec >= min(declim) & G3CGal$Dec <= max(declim) & G3CGal$Z >= min(zlim) & G3CGal$Z <= max(zlim))]),'CATAID'] 
		voidIDs = cat$voidgals[which(cat$voidgals[,'CATAID'] %in% G3CGal$CATAID[which(G3CGal$RA >= min(ralim) & G3CGal$RA <= max(ralim) & G3CGal$Dec >= min(declim) & G3CGal$Dec <= max(declim) & G3CGal$Z >= min(zlim) & G3CGal$Z <= max(zlim))]),'CATAID'] 

		#plot all groups in filaments
		
			plot(cat$groupcat[which(cat$groupcat[,'GroupID'] %in% filgroupIDs),'X'],cat$groupcat[which(cat$groupcat[,'GroupID'] %in% filgroupIDs),'Y'],cex=0.5,asp=1,axes=F,xlab="",ylab="",ylim=c(min(yrange)-10,max(yrange)+30),col='white',pch='.')
			if(pgals==T){points(cat$galcat[which(cat$galcat[,'CATAID'] %in% filgalIDs),'X'],cat$galcat[which(cat$galcat[,'CATAID'] %in% filgalIDs),'Y'],pch='.',col='blue',cex=4)}
			if(tendrils==T){points(cat$tendrilgals[which(cat$tendrilgals[,'CATAID'] %in% tendrilIDs),'X'],cat$tendrilgals[which(cat$tendrilgals[,'CATAID'] %in% tendrilIDs),'Y'],pch='.',col='green4',cex=4)}
			if(pvoids==T){points(cat$voidgals[which(cat$voidgals[,'CATAID'] %in% voidIDs),'X'],cat$voidgals[which(cat$voidgals[,'CATAID'] %in% voidIDs),'Y'],pch='.',col='red',cex=4)}


			if(mst==T){
				minimst = as.data.frame(cat$mst[which(cat$mst[,1] %in% filgroupIDs | cat$mst[,2] %in% filgroupIDs),])
				for(i in 1:dim(minimst)[1]){

					lines(x=c(cat$groupcat[which(cat$groupcat[,'GroupID']==as.numeric(as.character(minimst$From))[i]),'X'],cat$groupcat[which(cat$groupcat[,'GroupID']==as.numeric(as.character(minimst$To))[i]),'X']),y=c(cat$groupcat[which(cat$groupcat[,'GroupID']==as.numeric(as.character(minimst$From))[i]),'Y'],cat$groupcat[which(cat$groupcat[,'GroupID']==as.numeric(as.character(minimst$To))[i]),'Y']),col='cyan',lwd=1)

				}
			}
			overlaygrid(ralim=ralim,declim=declim,zlim=zlim,text=c(1,2,3),textsize=0.75,radeclabels=labels)
			
		}
}

#_________________________________________#
#_________________________________________#
