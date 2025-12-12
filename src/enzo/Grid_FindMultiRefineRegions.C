#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

void WriteListOfInts(FILE *fptr, int N, int nums[]);

int grid::FindMultiRefineRegions(int level)
{
    /* Return if this grid is not on this processor */
    if (MyProcessorNumber != ProcessorNumber)
        return SUCCESS;

    int region, i, CurrentMethod, j, CMind=INT_UNDEFINED;
    int UFind=0; // index of first LocalCellFlaggingMethod set to INT_UNDEFINED
    
    /* arrays to hold initial values for flagging methods and level limits */
    int IniLocalCellFlaggingMethod[MAX_FLAGGING_METHODS];
    int IniLocalMultiRefineMaximumLevel[MAX_FLAGGING_METHODS];
    int IniLocalMultiRefineMinimumLevel[MAX_FLAGGING_METHODS];

    /* Initialize local arrays to default values */
    for (i=0; i<MAX_FLAGGING_METHODS; i++){
        LocalCellFlaggingMethod[i] = INT_UNDEFINED;
        LocalMultiRefineMaximumLevel = MaximumRefinementLevel;
        LocalMultiRefineMinimumLevel = 0;
        IniLocalCellFlaggingMethod[i] = INT_UNDEFINED;
        IniLocalMultiRefineMaximumLevel[i] = INT_UNDEFINED;
        IniLocalMultiRefineMinimumLevel[i] = INT_UNDEFINED;
    }

    /* For each MultiRefineRegion, check whether there is overlap with this grid. If there is, 
       adjust the minimum and maximum refinement levels as necessary */
    for (region=0; region<NumberOfStaticMultiRefineRegions+NumberOfEnabledMultiRefineRegions; region++){
        /* If there is overlap between this grid and this MultiRefineRegion*/
        if(((GridRightEdge[0] > MultiRefineRegionLeftEdge[region][0]) && (GridLeftEdge[0] < MultiRefineRegionRightEdge[region][0]))
            && ((GridRightEdge[1] > MultiRefineRegionLeftEdge[region][1]) && (GridLeftEdge[1] < MultiRefineRegionRightEdge[region][1]))
            && ((GridRightEdge[2] > MultiRefineRegionLeftEdge[region][2]) && (GridLeftEdge[2] < MultiRefineRegionRightEdge[region][2]))){
            /* For each cell flagging method in this MultiRefineRegion...*/
            for (i = 0; i<MAX_FLAGGING_METHODS; i++){ 
                /* If this MultiRefineRegion is using this cell flagging slot...*/
                if (MultiRefineRegionFlaggingMethod[i] != INT_UNDEFINED){
                    CurrentMethod = MultiRefineRegionFlaggingMethod[i];
                    CMind = INT_UNDEFINED;
                    for (j=0; j<MAX_FLAGGING_METHODS; j++){
                        /* Identify the index of this cell flagging method (if it exists) */
                        if (IniLocalCellFlaggingMethod[j] == CurrentMethod){
                            CMind = j;
                        }
                    } // for (j=0; j<MAX_FLAGGING_METHODS; j++)
                    /* If this grid has not already been flagged for this method, flag it and adopt this region's
                       minimum and maximum refinement levels */
                    if (CMind==INT_UNDEFINED){ 
                        IniLocalCellFlaggingMethod[UFind] = CurrentMethod;
                        IniLocalMultiRefineMaximumLevel[UFind] = MultiRefineRegionMaximumLevel[region][i];
                        IniLocalMultiRefineMinimumLevel[UFind] = MultiRefineRegionMinimumLevel[region][i];
                        UFind++;
                    }
                    /* If this grid has already been flagged for this method, check what the minimum and maximum
                       permitted levels are for it and update them if necessary. The highest available value for
                       each will be used. */
                    else{
                        if (IniLocalMultiRefineMaximumLevel[CMind]<MultiRefineRegionMaximumLevel[region][i]){
                            IniLocalMultiRefineMaximumLevel[CMind] = MultiRefineRegionMaximumLevel[region][i];
                        }
                        if (IniLocalMultiRefineMinimumLevel[CMind]<MultiRefineRegionMinimumLevel[region][i]){
                            IniLocalMultiRefineMinimumLevel[CMind] = MultiRefineRegionMinimumLevel[region][i];
                        }
                    } // Does this grid already have this method?
                } // if (MultiRefineRegionFlaggingMethod[i] != INT_UNDEFINED)
            } // for (i = 0; i<MAX_FLAGGING_METHODS; i++)
        } // if this grid overlaps with this MultiRefineRegion
    } // for each MultiRefineRegion

    /* Check to see if there are any global cell flagging methods that are not yet flagged for this grid and
       add them. Minimum and maximum refinement levels adopt global values for a given refinement type or the
       general global values if more specific values have not been given */
    for (i = 0; i<MAX_FLAGGING_METHODS; i++){ 
        CMind = INT_UNDEFINED;
        for (j=0; j<MAX_FLAGGING_METHODS; j++){
            if (IniLocalCellFlaggingMethod[j] == CellFlaggingMethod[i]){
                CMind = j;
            }
        }
        if ((CMind == INT_UNDEFINED) && (CellFlaggingMethod[i] != INT_UNDEFINED)){
            IniLocalCellFlaggingMethod[UFind] = CellFlaggingMethod[i];
            if (CellFlaggingMethod[i] == 13){
                IniLocalMultiRefinementMinimumLevel[UFind] = MetallicityRefinementMinLevel;
            }
            else{
                IniLocalMultiRefineMinimumLevel[UFind] = 0;
            }
            if (CellFlaggingMethod[i] == 14){
                IniLocalMultiRefineMaximumLevel[UFind] = ShockwaveRefinementMaxLevel;
            }
            else{
                IniLocalMultiRefineMaximumLevel[UFind] = MaximumRefinementLevel;
            }
            UFind++;
        }
    }

    /* Find final MustRefine min and max for this region */
    for (i=0; i<MAX_FLAGGING_METHODS; i++){ 
        if (IniLocalCellFlaggingMethod[i] == 12){
            LocalMultiRefineMinimumLevel = IniLocalMultiRefineMinimumLevel[i];
            LocalMultiRefineMaximumLevel = IniLocalMultiRefineMaximumLevel[i];
        }
    }

    /* For each flagging method, check to see if the grid is already at or above the desired level.
       If it is, remove this flagging method. Note that if MustRefine is in use, the max specified
       by that method (i.e., LocalMultiRefineMaximumLevel) is the absolute maximum for this grid. */
    UFind = 0;
    for (i=0; i<MAX_FLAGGING_METHODS; i++){ 
        if (IniLocalCellFlaggingMethod[i] == 12){
            if (level < LocalMultiRefineMinimumLevel){
                LocalCellFlaggingMethod[UFind] = 12;
                UFind++;
            }
        }
        else{
            if ((level < IniLocalMultiRefineMaximumLevel[i]) && (level < LocalMultiRefineMaximumLevel)){
                LocalCellFlaggingMethod[UFind] = IniLocalCellFlaggingMethod[i];
                UFind++;

                /* if methods 13 and/or 14 are in use, update limit variables to local values */
                if (IniLocalCellFlaggingMethod[i] == 13){
                    MetallicityRefinementMinLevel = IniLocalMultiRefinementMinimumLevel[i];
                }
                if (IniLocalCellFlaggingMethod[i] == 14){
                    if (IniLocalMultiRefinementMaximumLevel[i] > LocalMultiRefineMaximumLevel){
                        ShockwaveRefinementMaxLevel = LocalMultiRefineMaximumLevel;
                    }
                    else{
                        ShockwaveRefinementMaxLevel = IniLocalMultiRefinementMaximumLevel[i];
                    }
                }
            }
        }

    }
    /* If no methods end up getting turned on, use method 0, which will set NumberOfFlaggedCells=0 for this grid */
    if (UFind==0){
        LocalCellFlaggingMethod[UFind] = 0;
    }

    if(debug1 && MyProcessorNumber == ROOT_PROCESSOR){
        fprintf(stderr, 'FindMultiRefineRegions says the following cell flagging methods have been turned on for this grid: ');
        WriteListOfInts(stderr, MAX_FLAGGING_METHODS, LocalCellFlaggingMethods);
    }
    

    return SUCCESS;
}