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
float *ini_met;
float *pop_age;

/*
    Arrays of size n_met * n_age * n_sources.
    n_sources is 4 for mass & metal fraction tables
    but 2 for SNe event rate table.
*/
float *mass_yield;
float *metm_yield;
float *event_rate;

};
