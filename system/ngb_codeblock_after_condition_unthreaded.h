if(P[p].Ti_current != ti_Current)
drift_particle(p, ti_Current);

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

Ngblist[numngb++] = p;
}
else
{
    if(no >= maxPart + maxNodes)	/* pseudo particle */
    {
        if(mode == 1) {endrun(123129);}
        if(target >= 0)	/* if no target is given, export will not occur */
        {
            if(Exportflag[task = DomainTask[no - (maxPart + maxNodes)]] != target)
            {
                Exportflag[task] = target;
                Exportnodecount[task] = NODELISTLENGTH;
            }
            
            if(Exportnodecount[task] == NODELISTLENGTH)
            {
                if(*nexport >= bunchSize)
                {
                    *nexport = nexport_save;
                    if(nexport_save == 0) {endrun(13004);} /* in this case, the buffer is too small to process even a single particle */
                    for(task = 0; task < NTask; task++) {nsend_local[task] = 0;}
                    for(no = 0; no < nexport_save; no++) {nsend_local[DataIndexTable[no].Task]++;}
                    return -1; /* buffer has filled -- important that only this and other buffer-full conditions return the negative condition for the routine */
                }
                Exportnodecount[task] = 0;
                Exportindex[task] = *nexport;
                DataIndexTable[*nexport].Task = task;
                DataIndexTable[*nexport].Index = target;
                DataIndexTable[*nexport].IndexGet = *nexport;
                *nexport = *nexport + 1;
                nsend_local[task]++;
            }
            
            DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]++] = DomainNodeIndex[no - (maxPart + maxNodes)];
            if(Exportnodecount[task] < NODELISTLENGTH) {DataNodeList[Exportindex[task]].NodeList[Exportnodecount[task]] = -1;}
        }
        no = Nextnode[no - maxNodes];
        continue;
    }
    
    current = &Nodes[no];
    
    if(mode == 1)
    {
        if(current->u.d.bitflags & (1 << BITFLAG_TOPLEVEL))	/* we reached a top-level node again, which means that we are done with the branch */
        {
            *startnode = -1;
#ifndef REDUCE_TREEWALK_BRANCHING
            return numngb;
#else
            return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml, SEARCHBOTHWAYS);
#endif
        }
    }
    
    if(current->Ti_current != ti_Current) {force_drift_node(no, ti_Current);}
    
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
    no = current->u.d.sibling;	// in case the node can be discarded //
#include "ngb_codeblock_checknode.h"
    no = current->u.d.nextnode;	// ok, we need to open the node //
}
}

*startnode = -1;
#ifndef REDUCE_TREEWALK_BRANCHING
return numngb;
#else
return ngb_filter_variables(numngb, Ngblist, &vcenter, &box, &hbox, hsml, SEARCHBOTHWAYS);
#endif

