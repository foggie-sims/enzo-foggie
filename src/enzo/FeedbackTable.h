/******************************************
/
/ Stellar Population Feedback Data Tables
/
/*****************************************/

struct FeedbackTableType {

/* 
    Each feedback table is 3 dimensions. 
    The first two are initial metal fraction and population age,
    while the third is the source.
    The arrays of metal fraction and population age are used to 
    interpolate these dimensions of the table.
    (See global_data.h for the source indexes.)
*/
int n_met;
int n_age;
double *ini_met;
double *pop_age;

/*
    Arrays of size n_met * n_age * n_sources.
    n_sources is 4 for mass & metal fraction tables
    but 2 for SNe event rate table.
*/
double *mass_yield;
double *metm_yield;
double *event_rate;

};

struct ChemFeedbackTableType {

/* 
    Chemical feedback table with the same dimensional structure
    as FeedbackTableType.
    The first two dimensions are initial metal fraction and population age,
    while the third is the element.
    The arrays of metal fraction and population age are used to 
    interpolate these dimensions of the table.
*/
int n_met;
int n_age;
double *ini_met;
double *pop_age;

/*
    Chemical yield arrays for individual elements.
    Each array is of size n_met * n_age * n_elem.
*/
double *C_yield;
double *O_yield;
double *Mg_yield;
double *Si_yield;
double *Fe_yield;

};
