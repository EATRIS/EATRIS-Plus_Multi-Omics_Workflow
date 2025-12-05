# check https://stackoverflow.com/questions/45385327/plotting-a-venn-digram-of-many-groups
## https://www.r-graph-gallery.com/14-venn-diagramm.html

# Venn = function(lst, main="", sub=""){
	# lst = lapply(lst,unique)
	# names(lst) = sprintf("%s: %d",names(lst),unlist(lapply(lst,length)))
	# v=venn(lst)
	# if (length(lst) == 2)main=sprintf("%s\nsimilarity = %.1f %%",main,100*2*v["11","num"]/sum(unlist(lapply(lst,length))))
	# title(main=main,sub=sub)
	# return(v)
# }
Venn = function(lst, main="", sub=""){
	require(gplots)
	names(lst) = sprintf("%s: %d",names(lst),unlist(lapply(lst,length)))
	lst = lapply(lst,unique)
	v=venn(lst)
	if (length(lst) == 2) jc = v["11","num"] / sum(v[,"num"])
	if (length(lst) == 3) jc = v["111","num"] / sum(v[,"num"])
	if (length(lst) == 4) jc = v["1111","num"] / sum(v[,"num"])
	
	main=sprintf("%s\nJaccard sim = %.4f",main,jc)
	title(main=main,sub=sub)
	return(v)
}

Euler = function(lst, main="", sub=""){
	require(VennDiagram)
	names(lst) = sprintf("%s: %d",names(lst),unlist(lapply(lst,length)))
	lst = lapply(lst,unique)
	
	if (length(lst) == 2){
		v = draw.pairwise.venn(area1 = length(lst[[1]]), 
							   area2 = length(lst[[2]]),
							   cross.area = sum(lst[[1]] %in% lst[[2]]),
							   category = names(lst),
							   rotation.degree=0) 
		jc = sum(lst[[1]] %in% lst[[2]])/(length(lst[[1]]) + length(lst[[1]]) - sum(lst[[1]] %in% lst[[2]]))
	}
	plot.new()
	grid.draw(v)
	title(main=sprintf("%s\nJaccard sim = %.4f",main,jc),sub=sub)
	return(v)
}


if (FALSE){
	Venn(list(A=letters[1:12],B=letters[5:20]))
	
	lst = list(A=letters[1:12],B=letters[5:20])
	main="Euler"
	sub="zzz"
	
	m = NULL
	for (i in 1:length(lst))
	unique(lst[1])
	
	v = venneuler(m)
	plot(venneuler(m))
	
}