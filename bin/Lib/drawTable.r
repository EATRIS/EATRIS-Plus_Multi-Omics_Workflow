################################################################################
## Draw a table to graphical window
##------------------------------------------------------------------------------
## data - data frame or matrix
## x0,y0 - starting coordinates
## dx,dy - vector or scalar width and height of cells
## row.names - set TRUE to draw row.names or FASLE to skip
## cex - font size
## col - vector or scalar colors of table cells
##------------------------------------------------------------------------------
## MIT license (c) P.Nazarov, 2011.  petr.nazarov[at]lih.lu
################################################################################

## ToDo: in new R problem appear for 1-column Data  !!!!!
## ToDo: add float formating

if (FALSE){ ## debug
data = iris
x0=0
y0=0.6
dx=c(0.1,0.1,0.1,0.1,0.1,0.25)
dy=0.05 
cex=1
row.names=TRUE
col=1
new=T
bg="white"
float.format = "%.2e"
}

drawTable = function(data,x0=0,y0=1,dx=0.2,dy=0.05,row.names=TRUE,cex=1, col=1, new=T, bg=NA, float.format = "%.2e",head.srt=0){
	if (class(data)[1]=="matrix") data = data.frame(data,check.names = F)
	## factor->character
	for (i in 1:ncol(data))
		if (class(data[,i])[1] == "factor") 
			data[,i]= as.character(data[,i])
	## Make character table Tab from data 
	Tab=data.frame(matrix(nr = nrow(data)+1, nc=ifelse(row.names,ncol(data)+1,ncol(data))))
	if (row.names) {
		Tab[-1,1] = rownames(data)
		Tab[1,1] = ""
		Tab[1,-1]=names(data)
		Tab[-1,-1]=data
		if (!is.null(float.format)){
			for (i in 1:ncol(data)){
				if (class(data[[i]])[1]%in%c("numeric","double")){
					if (class(Tab[-1,-1])[1] == "data.frame"){
						Tab[-1,-1][,i] = sprintf(float.format,data[[i]])
					}else{
						Tab[-1,-1] = sprintf(float.format,data[[i]])
					}
				}
			}
		}
	} else {
		Tab[1,]=names(data)
		Tab[-1,]=data
		if (!is.null(float.format))
			for (i in 1:ncol(data))
				if (class(data[[i]])[1]%in%c("numeric","double"))
					Tab[-1,][,i] = sprintf(float.format,data[[i]])
		## ToDo: add check for 1 column!
	}
	
	x = double(ncol(Tab))
	y = double(nrow(Tab))
	## check and elongate dx (if needed)
	if (length(dx) < length(x[-1])) 
		x = c(0,dx[-length(dx)],rep(dx[length(dx)],ncol(Tab)-length(dx)))
	if (length(dx) > length(x[-1])) 
		x = c(0,dx[1:(length(x)-1)])
	if (length(dx) == length(x[-1])) 
		x = c(0,dx)
	## check and elongate dy (if needed)
	if (length(dy) < nrow(Tab)-1) 
		y = c(0,dy[-length(dy)],rep(dy[length(dy)],nrow(Tab)-length(dy)))
	if (length(dy) > length(y[-1])) 
		y[-1] = dy[1: length(y[-1])]
	if (length(dy) == length(y[-1]))
		y[-1] = dy
	## check and elongate (in 2 dim) colors
	if (length(col) < ncol(Tab)*nrow(Tab)) {
		col = rep(col,ncol(Tab)*nrow(Tab))
	}
	
	x = cumsum(x) + x0
	y = y0 - cumsum(y)
	
	if (!new) 
	plot.new()
	#if (0) library(grDevices)
	k=1;iy=1;ix=1
	for (iy in 1:nrow(Tab)){
		if (y[iy]<mean(dy)){
			y[iy:nrow(Tab)] = y[iy:nrow(Tab)] + 1 - y[iy]
			plot.new()
		}
		for (ix in 1:ncol(Tab)){
			if (!is.na(bg)) rect(x[ix],y[iy], 1.05,y[iy]-mean(dy),col=bg,border=NA)
			text(x[ix],y[iy],
				Tab[iy,ix],
				adj=c(0,1),cex=cex,col=col[k],
				font = ifelse((ix==1 && row.names)|(iy==1),2,1),
				srt = ifelse(iy>1,0,head.srt))
			k=k+1
		}
		#getGraphicsEvent("click",onMouseDown = print("click"))
	}
}

# ## Example
# x11()
# plot.new()
# drawTable(USPersonalExpenditure,
          # dx=c(0.4,0.1,0.1,0.1,0.1),
          # dy = c(0.1,0.05,0.05,0.05,0.05),
          # row.names=T)
# # abline(h=c(0,1))
# # abline(v=c(0,1))
