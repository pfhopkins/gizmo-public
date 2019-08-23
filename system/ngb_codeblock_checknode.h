#if (BOX_SHEARING > 1)
/* in this case, we have a shearing box with the '1' coordinate being phi, so there is a periodic extra wrap; 
    this requires some extra care, because we can have the situation where the -node- center is wrapped, but 
    the actual particle[s] of interest are not wrapped */

    // start with the 'normal' periodic components: always pays to check these first because next is more expensive //
    //  (since there is no extra wrapping in these directions, the normal validation criteria apply)
    MyDouble dx0 = current->center[0]-searchcenter[0];
    dx = NGB_PERIODIC_BOX_LONG_X(dx0,current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
    if(dx > dist) continue;
    dz = NGB_PERIODIC_BOX_LONG_Z(dx0,current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
    if(dz > dist) continue;
    // ok we've made it this far, now we need to test the y-axis wrapping //
    MyDouble dx_abs = fabs(dx0);
    MyDouble dx_node = 0.5 * current->len;
    if((dx_abs + dx_node < boxHalf_X) || (dx_abs - dx_node > boxHalf_X))
    {
        // entire node box should be either wrapped or unwrapped: no ambiguity here, use the 'normal' distance evaluation
        dy = NGB_PERIODIC_BOX_LONG_Y(dx0,current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
        if(dy > dist) continue;
    } else {
        // ok this is the messy case: some of the node is wrapped, some is unwrapped. consider both cases
        MyDouble dy0 = current->center[1]-searchcenter[1];
        MyDouble dy_m = NGB_PERIODIC_BOX_LONG_Y(dx0-dx_node,dy0,current->center[2]-searchcenter[2],-1); // one edge of x-side of box
        MyDouble dy_p = NGB_PERIODIC_BOX_LONG_Y(dx0+dx_node,dy0,current->center[2]-searchcenter[2],-1); // other edge, should bracket possible wrapping
        dy = DMIN(dy_m, dy_p); // need to include all possible neighbors, which means using minimum of either distance
        if(dy > dist) continue;
    }
    // now test against the minimal sphere enclosing everything //
    dist += FACT1 * current->len;
    if(dx * dx + dy * dy + dz * dz > dist * dist) continue;

#else
/* this is the 'normal' operation mode */

#ifdef REDUCE_TREEWALK_BRANCHING
    // On the Power platform it is more efficient to compute all conditions first and then perform a single branch (on current Intel processors this is equally fast)
    MyDouble dist2 = dist + FACT1*current->len;
    dx = NGB_PERIODIC_BOX_LONG_X(current->center[0] - vcenter.d[0], current->center[1] - vcenter.d[1], current->center[2] - vcenter.d[2],-1);
    dy = NGB_PERIODIC_BOX_LONG_Y(current->center[0] - vcenter.d[0], current->center[1] - vcenter.d[1], current->center[2] - vcenter.d[2],-1);
    dz = NGB_PERIODIC_BOX_LONG_Z(current->center[0] - vcenter.d[0], current->center[1] - vcenter.d[1], current->center[2] - vcenter.d[2],-1);
    // now test against the minimal sphere enclosing everything
    if ((dx > dist) || (dy > dist) || (dz > dist) || (dx * dx + dy * dy + dz * dz > dist2 * dist2)) continue;    // ok, we need to open the node
#else
    // On older Intel and AMD processors it seems better to avoid the computations and branch early
    dx = NGB_PERIODIC_BOX_LONG_X(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
    if(dx > dist) continue;
    dy = NGB_PERIODIC_BOX_LONG_Y(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
    if(dy > dist) continue;
    dz = NGB_PERIODIC_BOX_LONG_Z(current->center[0]-searchcenter[0],current->center[1]-searchcenter[1],current->center[2]-searchcenter[2],-1);
    if(dz > dist) continue;
    // now test against the minimal sphere enclosing everything //
    dist += FACT1 * current->len;
    if(dx * dx + dy * dy + dz * dz > dist * dist) continue;
#endif
#endif
