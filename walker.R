require(fields)

# A bunch of functions designed to find backbones of filaments by walking around them from the edges to the centres. Written by Aaron Robotham, edited by Mehmet Alpaslan.

#_________________________________________#
#_________________________________________#

walk <- function(links){

    branchorder = 1
    links = as.matrix(links)
    cumwalk = cbind(unique(as.vector(links)), -1, 0)
    temp = table(c(links[, 1], links[, 2])) 
    # this table determines how many times each node is used.
    nodes = cbind(as.numeric(names(temp)), as.numeric(temp))
    Nfil = length(nodes[,1])
    start = nodes[nodes[, 2] == 1, 1] 
    # define starting points as nodes with only 1 edge.
    branchlist = list()
    newbranchlist = list()
    for (i in 1:length(start)) {
        count = 1
        previous = {}
        cumcount = {}
        branch = {}
        current = start[i]

        # Straightwalk is important. As long as it remains true, it means we are walking along a branch and haven't reached an intersection. The way to tell we've reached an intersection is for the table() count of the current command to be 1 (an edge) or > 2 (a branch).

        straightwalk = T
        while (straightwalk & count < Nfil/2) {
            cumwalk[cumwalk[, 1] == current, 2] = count # columns 2 and 3 of cumwalk defined as count and order
            cumwalk[cumwalk[, 1] == current, 3] = branchorder
            cumcount = c(cumcount, count)
            count = count + 1
            branch = c(branch, current)
            current = unique(as.vector(links[links[, 1] == current | links[, 2] == current, ])) 
            # add to current the next node along. either in from column or to column
            current = current[(!current %in% branch) & (!current %in% branch)] 
            # now search through previous and this branch for elements and remove them to change current node to the next one along.
            straightwalk = nodes[nodes[, 1] == current, 2] == 2 
            # straightwalk check to see if at edge, intersection, or keep going along branch
        }
        cumwalk[cumwalk[, 1] == current, 2] = count
        cumwalk[cumwalk[, 1] == current, 3] = branchorder # tally up again
        cumcount = c(cumcount, count)
        branch = c(branch, current)

        # Everything goes into branchlist. Fill it up as you go along (goodbye Fortran allocatable arrays!)

        branchlist[[i]] = list(branch = branch, count = count, cumcount = cumcount, previous = previous) 
    }

    # Now all edges have been made into branches until reaching an intersection. Next step is to look at all nodes with count = -1 (i.e. not visited yet). If, however, all points have already been visited, then set the end of the longest branch as the centre and output.

    # This can be a problem. If all nodes are visited and there are only 2 branches, then just output. However, in some cases all nodes can be visited but there are many branches (complex filaments) and in those cases you must continue.

    if(all(cumwalk[,2]>-1) & length(branchlist) == 2){
        branchlengths = {}
        Nbranch = length(branchlist)
        for(i in 1:Nbranch){
            branchlengths = append(branchlengths,branchlist[[i]]$count)
        }
        
        cumwalk[cumwalk[, 1] == rev(branchlist[[which.max(branchlengths)]]$branch)[1], 2] = length(cumwalk[, 1])

        return = cumwalk = cumwalk

    }else{

    ends = c(1,2,3);tomerge=1
    while (length(which(cumwalk[, 2] == -1)) > 0 | length(unique(ends)) > 2) {
        branchorder = branchorder + 1
        ends = {}
        Nbranch = length(branchlist)

        # Start by merging branches with a common intersection.

        for (i in 1:Nbranch) {
            ends = c(ends, rev(branchlist[[i]]$branch)[1]) 
            # ends are found by reversing branches.
        }
        temp = table(ends)
        tomerge = cbind(as.numeric(names(temp)), as.numeric(temp))
        tomerge = tomerge[tomerge[, 2] > 1, 1]  
        # again, check table() entry to find which endpoints share multiple branches and choose to merge those. 
        merged = {}
        combsum = {}

        # Go through each end that appears more than once and sum branch elements.

        for (j in 1:length(tomerge)) {
            mergebranch = which(ends == tomerge[j])
            merged = c(merged, mergebranch) # save IDs of merged branches here.
            sumsofar = 0
            for (k in 1:length(mergebranch)) {
                sumsofar = sumsofar + branchlist[[mergebranch[k]]]$count
            }
            combsum = c(combsum, sumsofar)
        }
        keep = 1:Nbranch
        keep = keep[!keep %in% merged] 
        # keep all branches that have not been merged. Of all branches (1:Nbranch), which ones aren't in the list of merged branches?
        tomerge = tomerge[order(combsum)]

        # Now repeat the process of finding next node along branches until you reach another interesction. 

        for (j in 1:length(tomerge)) {
            mergebranch = which(ends == tomerge[j])
            count = 1
            previous = {}

            # Tally up count from precursor branches.

            for (k in 1:length(mergebranch)) {
                previous = c(previous, branchlist[[mergebranch[k]]]$branch, branchlist[[mergebranch[k]]]$previous)
            }

            # Keep track of previous branches for being able to move current along.

            previous = unique(previous)
            count = length(previous) - 1
            cumcount = {}
            branch = {}
            current = tomerge[j]

            # Now do the straightwalk branch progression again, as before.

            straightwalk = T
            while (straightwalk & count < Nfil/2) {
                cumwalk[cumwalk[, 1] == current, 2] = count
                cumwalk[cumwalk[, 1] == current, 3] = branchorder
                cumcount = c(cumcount, count)
                count = count + 1
                branch = c(branch, current)
                current = unique(as.vector(links[links[, 1] == current | links[, 2] == current, ]))
                current = current[(!current %in% previous) & (!current %in% branch)]
                straightwalk = length(current) == 1
            }
            if (length(current) == 1) {
                cumcount = c(cumcount, count)
                branch = c(branch, current)
            }

            # Make a new branch list for this next phase.

            newbranchlist[[j]] = list(branch = branch, count = count, cumcount = cumcount, previous = previous)
        }
        if (length(keep) > 0) {
            for (i in 1:length(keep)) {
                newbranchlist = c(newbranchlist, list(branchlist[[keep[i]]]))
            }
        }
        branchlist = newbranchlist
        newbranchlist = list()

    }
    # Keep going until all nodes visited, then do one last check for objects to merge.

    ends = {}
    Nbranch = length(branchlist)
    for (i in 1:Nbranch) {
            ends = c(ends, rev(branchlist[[i]]$branch)[1]) 
        }
    temp = table(ends)
    tomerge = cbind(as.numeric(names(temp)), as.numeric(temp))
    tomerge = tomerge[tomerge[, 2] > 1, 1]
    merged = {}
    combsum = {}

    if(length(tomerge) > 0){ #only do this if tomerge has elements in it, otherwise will crash.

        for (j in 1:length(tomerge)) {
            mergebranch = which(ends == tomerge[j])
            merged = c(merged, mergebranch) # save IDs of merged branches here.
            sumsofar = 0
            for (k in 1:length(mergebranch)) {
                sumsofar = sumsofar + branchlist[[mergebranch[k]]]$count
            }
            combsum = c(combsum, sumsofar)
        }

        # Final check for objects to merge.

        for (j in 1:length(tomerge)) {
                mergebranch = which(ends == tomerge[j])
                count = 1
                previous = {}

                # Tally up count from precursor branches.

                for (k in 1:length(mergebranch)) {
                    previous = c(previous, branchlist[[mergebranch[k]]]$branch, branchlist[[mergebranch[k]]]$previous)
                }

                # Keep track of previous branches for being able to move current along.

                previous = unique(previous)
                count = length(previous) - 1
                cumcount = {}
                branch = {}
                current = tomerge[j]
                cumwalk[cumwalk[, 1] == current, 2] = count
                cumwalk[cumwalk[, 1] == current, 3] = branchorder
        }
    }

    # Set count for central node to the max, then output. Central node will be at the end of the longest filament.

    branchlengths = {}
    Nbranch = length(branchlist)
    for(i in 1:Nbranch){
        branchlengths = append(branchlengths,branchlist[[i]]$count)
    }
    
    cumwalk[cumwalk[, 1] == rev(branchlist[[which.max(branchlengths)]]$branch)[1], 2] = length(cumwalk[, 1])

    return = cumwalk = cumwalk}
}


#_________________________________________#
#_________________________________________#

makefilament <- function (inwalk, links, start = "max", doavoid = T, direction = "down", dodist = F, xyz = {}){

    if (dodist & length(xyz) == 0) {
        stop("Need xyz!!")
    }

    # Starting at max means just starting at the node which has the largest count value from walk()'s output. i.e. node that is furthest from the edges.

    maxloc = inwalk[which.max(inwalk[, 2]), 1]
    if (start[1] == "max") {
        start = maxloc
        direction = "down"
    }

    if (length(start) == 1) {
        doavoid = T

        # posdir identifies all nodes that are linked to the starting point.

        posdir = unique(as.vector(links[links[, 1] %in% start | links[, 2] %in% start, ]))
        posdir = posdir[!posdir %in% start]
        Nroutes = length(posdir)

        # routes then goes through each of these.

        routes = list()
        avoid = {}

        for (i in 1:Nroutes) {
            filament = start
            position = start
            cont = T
            if (direction == "down") {
                filcount = 1e+10
            }
            if (direction == "up") {
                filcount = -1e+10
            }
            # cont is the flag that determines if the finder continues or stops.
            while (cont == T) {
                # Set up posdir such that it avoids objects in filament, and anything to avoid.

                posdir = unique(as.vector(links[links[, 1] %in% position | links[, 2] %in% position, ]))
                posdir = posdir[(!posdir %in% filament) & (!posdir %in% avoid)]
                if (length(posdir) > 0) {
                  if (direction == "down") {

                    # check the walk() output for each neighbouring node and move onto the one with the next highest count, such that the longest filament is made.
                    if (max(inwalk[match(posdir, inwalk[, 1]), 2]) <= filcount) {
                      position = posdir[which.max(inwalk[match(posdir, inwalk[, 1]), 2])]
                      filcount = max(inwalk[match(posdir, inwalk[,1]), 2])
                      if (doavoid) {
                        # if doavoid is true, the avoid vector is filled in with the current position. This ensures it cannot be revisited.
                        avoid = c(avoid, position)
                      }
                      filament = c(filament, position)

                      # If the filament has reached the centre, or an edge, stop.
                      if (inwalk[inwalk[, 1] == position, 2] == 0 | position == maxloc) {
                        cont = F
                      }
                    }
                    else {
                      cont = F
                    }
                  }
                  if (direction == "up") {

                    # Same as direction = down, except aim for the node with the highest count.
                    if (max(inwalk[match(posdir, inwalk[, 1]), 2]) >= filcount) {
                      position = posdir[which.max(inwalk[match(posdir, inwalk[, 1]), 2])]
                      filcount = max(inwalk[match(posdir, inwalk[, 1]), 2])
                      if (doavoid) {
                        avoid = c(avoid, position)
                      }
                      filament = c(filament, position)
                      if (inwalk[inwalk[, 1] == position, 2] == 0 | position == maxloc) {
                        cont = F
                      }
                    }
                    else {
                      cont = F
                    }
                  }
                }
                else {

                # This is important. The branch has been constructed but it does not link to a parent. This is not a big deal but it's important for distances to be calculated. So set up a temporary filament to put into linkdist below, which contains the branch connection to the parent.

                     if(position != maxloc){
                        posdir = unique(as.vector(links[links[, 1] %in% position | links[, 2] %in% position, ]))
                        posdir = posdir[(!posdir %in% filament)]
                        filament_temp={};filament_temp = c(filament,posdir)
                        cont = F
                    }
                  
                }
            }

            if(exists("filament_temp")==T){linksback = rbind(links[links[, 1] %in% filament_temp & links[, 2] %in% filament_temp, ])}else{linksback = rbind(links[links[, 1] %in% filament & links[, 2] %in% filament, ])}

            # if dodist = T, this chunk of code computes distances between links and adds them to the links data frame.
            if (dodist & length(linksback) > 0) {
                linkdist = {
                }
                linkxyz = cbind(xyz[match(linksback[, 1], xyz[, 1]), 2:4], xyz[match(linksback[, 2], xyz[,1]), 2:4])
                for (j in 1:length(linkxyz[, 1])) {
                  linkdist = c(linkdist, rdist(linkxyz[j, 1:3], linkxyz[j, 4:6]))
                }
                linksback = cbind(linksback, linkdist)
            }
            else {
                linksback = cbind(NA, NA, NA)
            }

            # Add this filament to the routes list.
            routes[[i]] = list(fil = filament, links = linksback, N = length(filament))
        }
    }

    # The above is for filaments that start from one place. Now to look at what to do when starting from multiple places. 

    # Mostly the same code, except the number of routes is given by the number of starting points, and you go through them all one by one.

    else {
        Nroutes = length(start)
        posdir = unique(as.vector(links[links[, 1] %in% start | links[, 2] %in% start, ]))
        posdir = posdir[!posdir %in% start]
        routes = list()
        avoid = {
        }
        for (i in 1:Nroutes) {
            filament = start[i]
            position = start[i]
            cont = T
            if (direction == "down") {
                filcount = 1e+10
            }
            if (direction == "up") {
                filcount = -1e+10
            }
            while (cont == T) {
                posdir = unique(as.vector(links[links[, 1] %in% position | links[, 2] %in% position, ]))
                posdir = posdir[(!posdir %in% filament) & (!posdir %in% avoid)]
                if (length(posdir) > 0) {
                  if (direction == "down") {
                    if (max(inwalk[match(posdir, inwalk[, 1]), 2]) <= filcount) {
                      position = posdir[which.max(inwalk[match(posdir, inwalk[, 1]), 2])]
                      filcount = max(inwalk[match(posdir, inwalk[,1]), 2])
                      if (doavoid) {
                        avoid = c(avoid, position)
                      }
                      filament = c(filament, position)
                      if (inwalk[inwalk[, 1] == position, 2] == 0 | position == maxloc) {
                        cont = F
                      }
                    }
                    else {
                      cont = F
                    }
                  }
                  if (direction == "up") {
                    if (max(inwalk[match(posdir, inwalk[, 1]), 2]) >= filcount) {
                      position = posdir[which.max(inwalk[match(posdir, inwalk[, 1]), 2])]
                      filcount = max(inwalk[match(posdir, inwalk[, 1]), 2])
                      if (doavoid) {
                        avoid = c(avoid, position)
                      }
                      filament = c(filament, position)
                      if (inwalk[inwalk[, 1] == position, 2] == 0 | position == maxloc) {
                        cont = F
                      }
                    }
                    else {
                      cont = F
                    }
                  }
                }
                else {
                    if(position != maxloc){
                        posdir = unique(as.vector(links[links[, 1] %in% position | links[, 2] %in% position, ]))
                        posdir = posdir[(!posdir %in% filament)]
                        filament_temp={};filament_temp = c(filament,posdir)
                        cont = F
                    }
                }
            }

            if(exists("filament_temp")==T){linksback = rbind(links[links[, 1] %in% filament_temp & links[, 2] %in% filament_temp, ])}else{linksback = rbind(links[links[, 1] %in% filament & links[, 2] %in% filament, ])}
            if (dodist & length(linksback) > 0) {
                linkdist = {}
                linkxyz = cbind(xyz[match(linksback[, 1], xyz[, 1]), 2:4], xyz[match(linksback[, 2], xyz[,1]), 2:4])
                for (j in 1:length(linkxyz[, 1])) {
                  linkdist = c(linkdist, rdist(linkxyz[j, 1:3], linkxyz[j, 4:6]))
                }
                linksback = cbind(linksback, linkdist)
            }
            else {
                linksback = cbind(NA, NA, NA)
            }
            routes[[i]] = list(fil = filament, links = linksback, N = length(filament))
        }
    }
    return = routes
}

#_________________________________________#
#_________________________________________#

findback <- function (links, bydist = F, dodist = F, xyz = {}){

    if (bydist == T) {
        dodist = T
    }

    # Start by running walk and makefilament -- starting it from edges, save outputs. Check records where edges are in walk.

    temp = walk(links)
    check = temp[temp[, 2] == 1, 1]
    fils = makefilament(temp, links, check, doavoid = F, direction = "up", dodist = dodist, xyz = xyz)
    Nfil = length(fils)
    fillen = {}

    # Record sizes of branches.

    for (i in 1:Nfil) {
        fillen = c(fillen, fils[[i]]$N)
    }

    # If xyz distance doesn't matter, reorder check to go from longest to shortest filament in terms of nodes. Else, go by cartesian distance.

    if (bydist == F) {
        check = check[order(fillen, decreasing = T)]
    }
    else {
        fildist = {}
        for (i in 1:Nfil) {
            fildist = c(fildist, sum(fils[[i]]$links[, 3]))
        }
        check = check[order(fildist, decreasing = T)]
    }

    # Rerun makefilament with ordered checks, this time avoiding. Since we know where the biggest branches start, by rerunning makefilament with doavoid = T and the starting positions going in descending order, we ensure the top branches will meet in the middle.

    fils = makefilament(temp, links, check, doavoid = T, direction = "up", dodist = dodist, xyz = xyz)
    Nfil = length(fils)
    fillen = {}

    for (i in 1:Nfil) {
        fillen = c(fillen, fils[[i]]$N)
    }

    # Again, if interested in distances, calculate cartesian distances.

    if (dodist) {
        fildist = {}
        for (i in 1:Nfil) {
            fildist = c(fildist, sum(fils[[i]]$links[, 3]))
        }
    }

    # Now go through and populate temp2, a list, with branches in decreasing size of number of nodes.

    temp2 = list()
    if (bydist == F) {
        for (i in 1:Nfil) {
            temp2[[i]] = fils[[order(fillen, decreasing = T)[i]]]
        }
        if (dodist) {
            fildist = fildist[order(fillen, decreasing = T)]
        }
        fillen = sort(fillen, decreasing = T)
    }

    # If doing by distance, then arrange by longest distance branches.

    else {
        for (i in 1:Nfil) {
            temp2[[i]] = fils[[order(fildist, decreasing = T)[i]]]
        }
        fillen = fillen[order(fildist, decreasing = T)]
        fildist = sort(fildist, decreasing = T)
    }

    # Then output! Only difference here is that if bydist = T, output list has an additional element of distances.

    if (dodist == F) {
        return = list(fils = temp2, Nfil = fillen, walk = temp)
    }
    else {
        return = list(fils = temp2, Nfil = fillen, walk = temp, dist = fildist)
    }
}

makenetwork <- function(links,bydist=F,dodist=F,xyz={}){

    # This function takes the output of findback and constructs the filamentary network by identifying the backbone and its tributaries.

    # Throw an rbind() around links to make sure it works even if there are 2 groups in the filament.

    links = rbind(links)

    fb = findback(links,bydist,dodist,xyz)

    # Now construct backbone and network.
    # Start by identifying centre, ends, and neighbouring nodes to centre.

    centre = rev(fb$fils[[1]]$fil)[1]
    midlinks = links[c(which(links[,1]==centre),which(links[,2]==centre)),]
    neighbours = midlinks[which(midlinks!=centre)]

    # Pick out second longest filament which joins to the centre at one end and has a endpoint at the other.

    i = 2
    while(any(fb$fils[[i]]$fil %in% neighbours)==F){i=i+1}
    # i is now set to the filament which fulfills this criterion.

    backbone = c(fb$fils[[1]]$fil,rev(fb$fils[[i]]$fil))
    Nbb = length(backbone)

    # save i for later. So filcounter is the id of the branch that joins to branch 1 to form the backbone. Toss this into branchlist later.

    filcounter = i

    # Make links.
    bblinks = {}

    for(i in 1:(Nbb-1)){

        bblinks = rbind(bblinks,cbind(backbone[i],backbone[i+1]))

    }

    # Calculate total distance of backbone.

    if(dodist & length(bblinks) > 0){
    
        bbdist = {}
        bbxyz = cbind(xyz[match(bblinks[, 1], xyz[, 1]), 2:4], xyz[match(bblinks[, 2], xyz[,1]), 2:4])
        for (j in 1:length(bbxyz[, 1])) {
            bbdist = c(bbdist, rdist(bbxyz[j, 1:3], bbxyz[j, 4:6]))
        }

        bblinks = cbind(bblinks,as.vector(bbdist))
        bblist = list(backbone = backbone, Nbb = Nbb, bblinks = bblinks, bbdist = sum(bbdist))
    }else{
        bblist = list(backbone = backbone, Nbb = Nbb, bblinks = bblinks)
    }

    # Now go through the remaining branches and assign orders to them, based on if they are neighbours of the backbone, or neighbours of sub-branches, etc...

    newmidlinks = links[links[,1] %in% backbone | links[,2] %in% backbone,]
    newneighbours = newmidlinks[which(newmidlinks %in% backbone==F)]
    filorder = 2

    # Loop over branches until no new neighbours left. After each set of branches, increase order.
    branchlist = {}

    # Add the branches that form the backbone as order number 1 branches.

    branchlist = rbind(branchlist,c(1,1))
    branchlist = rbind(branchlist,c(filcounter,1))

    # Continue!
    
    newnewneighbours = {}
    previous = {}
    while(length(newneighbours) > 0){

        for(j in 1:length(newneighbours)){
            i = 2
            while(any(fb$fils[[i]]$fil %in% newneighbours[j])==F){i=i+1}

            # New branch identified. Add this to branchlist and populate temp3, which will store new new neighbours. branchlist will be appended with the numbers of branches that connect to the backbone. From this identify new neighbours and continue until no neighbours are found.
            branchlist = rbind(branchlist,c(i,filorder))
            previous = c(previous,fb$fils[[i]]$fil)

            # Temp3 will contain new neighbours to the branch.
            temp3 = links[links[,1] %in% fb$fils[[i]]$fil | links[,2] %in% fb$fils[[i]]$fil,]
            temp3 = temp3[temp3 %in% fb$fils[[i]]$fil==F]
            temp3 = temp3[temp3 %in% previous==F]
            temp3 = temp3[temp3 %in% backbone==F]
            if(length(temp3)==0){next}
            newnewneighbours = c(newnewneighbours,as.numeric(temp3))
        }

        newneighbours = newnewneighbours
        newnewneighbours = {}
        filorder = filorder + 1
    }

    if(length(branchlist)==0){branchlist=matrix(c(1,1),nrow=1,ncol=2)}
    # Output.

    colnames(branchlist) = c("BranchID","BranchOrder")

    outlist = list(bblist=bblist, network = branchlist)
    return(c(fb,outlist))

}
