/******************************************
/
/ Stellar Population Feedback Data Tables
/
/*****************************************/

struct PreSNFeedbackTableType {

/* 
    Each pre-SN feedback table is 2 dimensions. 
    The dimensions are initial metal fraction and population age.
    The arrays of metal fraction and population age are used to 
    interpolate these dimensions of the table.
*/
int n_met;
int n_age;
double *ini_met;
double *pop_age;

/*
    Arrays of size n_met * n_age.
*/
double *mass_yield;
double *metm_yield;
double *mom_rate;

};
