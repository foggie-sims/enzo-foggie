/***********************************************************************
/
/  PROTO SUBGRID CLASS
/
/  written by: Greg Bryan
/  date:       October, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#ifndef PROTO_SUBGRID_DEFINED__
#define PROTO_SUBGRID_DEFINED__

class ProtoSubgrid
{
 private:

  int GridRank;
  int GridDimension[MAX_DIMENSION];
  
  int Level;

  FLOAT GridLeftEdge[MAX_DIMENSION];
  FLOAT GridRightEdge[MAX_DIMENSION];

  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  int NumberFlagged;

  int  *GridFlaggingField;
  int  *Signature[MAX_DIMENSION];

  /* Work-capped Berger-Rigoutsos (SubgridMaximumWorkFraction).  A
     ProtoSubgrid already indexes into the parent's arrays via
     StartIndex/GridDimension, so its predicted work can be summed
     straight from the parent's WorkField - no parallel array has to be
     copied through every split.  Both are NULL when the feature is off. */
  float *ParentWorkField;
  int    ParentDimension[MAX_DIMENSION];

 public:

  ProtoSubgrid();
  ~ProtoSubgrid();

  int AcceptableSubgrid();
  int ReturnNthLongestDimension(int n);
  int ComputeSignature(int dim);
  int FindGridsByZeroSignature(int dim, int &NumberOfNewGrids, 
			       int GridEnds[MAX_NUMBER_OF_SUBGRIDS][2]);
  int CopyToNewSubgrid(int dim, int GridStart, int GridEnd, 
		       ProtoSubgrid *NewGrid);
  int ComputeSecondDerivative(int dim, int &ZeroCrossStrength, 
			      int GridEnds[2][2]);
  int LargeAxisRatioCheck(int &dim, int GridEnds[MAX_DIMENSION*2][2], 
			  float CriticalRatio);
  int CopyFlaggedZonesFromGrid(grid *Grid);
  float ComputeWork();
  int ShrinkToMinimumSize();
  int CleanUp();

  int ReturnGridRank() {return GridRank;};
  int *ReturnGridDimension() {return GridDimension;};
  FLOAT *ReturnGridLeftEdge() {return GridLeftEdge;};
  FLOAT *ReturnGridRightEdge() {return GridRightEdge;};
  
  /* Return, set level */

  void SetLevel(int level) { Level = level; };
  int GetLevel(void) { return Level; };

};


#endif
