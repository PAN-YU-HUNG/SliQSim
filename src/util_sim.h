#ifndef _UTIL_SIM_H_
#define _UTIL_SIM_H_

#include <iostream>

/* function */
extern void full_adder_plus_1(int length, int *reg);
extern void full_adder_plus_1_start(int length, int *reg, int start);
extern void full_adder_plus_1_measure(int length, int *reg, int *order);
extern void replace_all(std::string &str, const std::string &search_str, const std::string &replace_str);
extern int int_array_full_check(int length, int *reg);
extern size_t getPeakRSS();
extern size_t getCurrentRSS();
extern std::string find_in_yaml(std::istream& yaml, const std::string& key);
extern double parse_theta(std::string theta);
extern double normalize_theta(double theta);

#endif
