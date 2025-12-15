/* Struct for variables relating to individual time entries 
   for individual MultiRefineRegion Tracks; member of 
   MultiRefineRegionTrack struct (see MultiRefineRegionTrack.h).
   Variables are populated in ReadEvolveRefineFile.C
*/

struct MultiRefineRegionTimeEntry {

    /* Time for this time entry
    If MultiRefineRegionTimeType == 0, this is code time 
    If MultiRefineRegionTimeType == 1, this is redshift */
    double Time;

    /* Holds coordinates of lower left corner and upper right
    corner of MultiRefineRegion at this time */
    double Pos[6]; 

    /* Arrays that contain list of minimum and maximum levels
    permitted for MultiRefineRegion at this time. Order is
    the same as MultiRefineRegionTrack.RefTypes */
    int *MinLevels;
    int *MaxLevels;

    /* Minimum star particle mass that can be formed in this
    MultiRefineRegion at this time in Msol */
    double MinStarMass;
};