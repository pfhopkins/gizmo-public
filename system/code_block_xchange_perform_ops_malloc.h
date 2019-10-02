/*! allocate buffers to arrange communication */
long long NTaskTimesNumPart = maxThreads * NumPart; size_t MyBufferSize = All.BufferSize; int loop_iteration = 0;
All.BunchSize = (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) + sizeof(struct INPUT_STRUCT_NAME) + sizeof(struct OUTPUT_STRUCT_NAME) + sizemax(sizeof(struct INPUT_STRUCT_NAME),sizeof(struct OUTPUT_STRUCT_NAME))));
Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));
DataIndexTable = (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
DataNodeList = (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));
double timeall=0, timecomp=0, timecomm=0, timewait=0, t0; CPU_Step[CPU_MISC] += measure_time(); t0 = my_second(); /*! for timing information */

