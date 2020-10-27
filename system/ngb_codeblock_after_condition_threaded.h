/*
 this defines a code-block to be inserted in the neighbor search routines after the conditions for neighbor-validity are applied
 (valid particle types checked)
 */
if(P[p].Ti_current != ti_Current)
{
    LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
    drift_particle(p, ti_Current);
    UNLOCK_PARTNODEDRIFT;
}

#ifndef REDUCE_TREEWALK_BRANCHING
#if (SEARCHBOTHWAYS==1)
dist = DMAX(PPP[p].Hsml, hsml);
#else
dist = hsml;
#endif
dx = NGB_PERIODIC_BOX_LONG_X(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
if(dx > dist) continue;
dy = NGB_PERIODIC_BOX_LONG_Y(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
if(dy > dist) continue;
dz = NGB_PERIODIC_BOX_LONG_Z(P[p].Pos[0] - searchcenter[0], P[p].Pos[1] - searchcenter[1], P[p].Pos[2] - searchcenter[2],-1);
if(dz > dist) continue;
if(dx * dx + dy * dy + dz * dz > dist * dist) continue;
#endif
ngblist[numngb++] = p;  /* Note: unlike in previous versions of the code, the buffer can hold up to all particles. note also the threaded-vs-unthreaded use of n vs N in ngblist */
}
else
{
    if(no >= maxPart + maxNodes)	/* pseudo particle */
    {
#ifdef DONOTUSENODELIST
        if(mode == 1)
        {
            no = Nextnode[no - maxNodes];
            continue;
        }
#endif
        if(mode == 1) {endrun(123128);}
        
        if(target >= 0)	/* if no target is given, export will not occur */
        {
            if(exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
            {
                exportflag[task] = target;
                exportnodecount[task] = NODELISTLENGTH;
            }
            
            if(exportnodecount[task] == NODELISTLENGTH)
            {
                int exitFlag = 0, nexp;
                LOCK_NEXPORT;
#ifdef _OPENMP
#pragma omp critical(_nexport_)
#endif
                {
                    if(Nexport >= bunchSize)
                    {
                        /* out of buffer space. Need to discard work for this particle and interrupt */
                        BufferFullFlag = 1;
                        exitFlag = 1;
                    }
                    else
                    {
                        nexp = Nexport;
                        Nexport++;
                    }
                }
                UNLOCK_NEXPORT;
                if(exitFlag) {return -1;} /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
                
                exportnodecount[task] = 0;
                exportindex[task] = nexp;
                DataIndexTable[nexp].Task = task;
                DataIndexTable[nexp].Index = target;
                DataIndexTable[nexp].IndexGet = nexp;
            }
#ifndef DONOTUSENODELIST
            DataNodeList[exportindex[task]].NodeList[exportnodecount[task]++] = DomainNodeIndex[no - (maxPart + maxNodes)];
            if(exportnodecount[task] < NODELISTLENGTH)
                DataNodeList[exportindex[task]].NodeList[exportnodecount[task]] = -1;
#endif
                }
        
        no = Nextnode[no - maxNodes];
        continue;
    }
    
    current = &Nodes[no];
    
#ifndef DONOTUSENODELIST
    if(mode == 1)
    {
        if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
        {
            *startnode = -1;
#ifndef REDUCE_TREEWALK_BRANCHING
            return numngb;
#else
            return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml, SEARCHBOTHWAYS);
#endif
        }
    }
#endif
    
    if(current->Ti_current != ti_Current)
    {
        LOCK_PARTNODEDRIFT;
#ifdef _OPENMP
#pragma omp critical(_partnodedrift_)
#endif
        force_drift_node(no, ti_Current);
        UNLOCK_PARTNODEDRIFT;
    }
    
    if(!(current->u.d.bitflags & (1 << BITFLAG_MULTIPLEPARTICLES)))
    {
        if(current->u.d.mass)	/* open cell */
        {
            no = current->u.d.nextnode;
            continue;
        }
    }
    
#if (SEARCHBOTHWAYS==1)
    dist = DMAX(Extnodes[no].hmax, hsml) + 0.5 * current->len;
#else
    dist = hsml + 0.5 * current->len;
#endif
    no = current->u.d.sibling;	/* in case the node can be discarded */
#include "ngb_codeblock_checknode.h"
    no = current->u.d.nextnode;	// ok, we need to open the node //
}
}

*startnode = -1;
#ifndef REDUCE_TREEWALK_BRANCHING
return numngb;
#else
return ngb_filter_variables(numngb, ngblist, &vcenter, &box, &hbox, hsml, SEARCHBOTHWAYS);
#endif
