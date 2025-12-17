/* Struct for variables relating to individual MultiRefineRegion 
   tracks. Variables are populated in ReadEvolveRefineFile.C */

#include "MultiRefineRegionTimeEntry.h"

struct MultiRefineRegionTrackType {

    /* Is this track currently enabled? 0 for no, 1 for yes */
    int Enabled;

    /* Number of flagging methods for which minimum and maximum
       levels are specified for this MultiRefineRegion */
    int NRefTypes; 

    /* List of flagging methods specified for this MultiRefineRegion */
    int *RefTypes;

    /* Number of time entries over which positions, min/max levels, etc
       are specified for this MultiRefineRegion */
    int NTimeEntries;

    /* Struct for variables related to individual time entries for this 
       MultiRefineRegion - see MultiRefineRegionTimeEntry.h */
    MultiRefineRegionTimeEntryType *TimeEntries;
};