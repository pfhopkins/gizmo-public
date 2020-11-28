/****************************************************************************
 * This file is part of CHIMES.
 * Copyright (c) 2020 Alexander Richings (alexander.j.richings@durham.ac.uk)
 *
 * CHIMES is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ***************************************************************************/

#ifdef CHIMES_ENABLE_GNU_SOURCE
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 
#endif 
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <time.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h> 
#include <sundials/sundials_types.h>
#include <sunlinsol/sunlinsol_dense.h> 
#include <sunmatrix/sunmatrix_dense.h> 
#include "../../allvars.h"
#include "chimes_vars.h"


