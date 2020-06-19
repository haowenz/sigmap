#ifndef FAST_DTW_H_
#define FAST_DTW_H_

#include<iostream>
#include<fstream>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include<string>

#include "signal_batch.h"
#include "utils.h"

namespace sigmap {

struct Coordinate{
   size_t target_coord;
   size_t query_coord;
   Coordinate(size_t new_target_coord,size_t new_query_coord):
    target_coord(new_target_coord),query_coord(new_query_coord){}
    bool operator < (const Coordinate &other) const { return ((target_coord < other.target_coord)||((target_coord == other.target_coord)&(query_coord < other.query_coord))); }
};
// std::ostream& operator<<(std::ostream& os, const Coordinate& c){
//     os <<  c.target_coord << "," << c.query_coord; 
//     return os;
// }
void reduce_by_half(const float *signal, size_t length, float *signal_reduced, size_t &reduced_length);
void expand_window(std::vector<Coordinate> &path,size_t target_length,size_t query_length,int radius,std::vector<std::vector<Coordinate>> &window);
void generate_path(std::vector<std::vector<int>> &path_matrix,std::vector<std::vector<Coordinate>> &window,std::map<Coordinate,std::pair<size_t,size_t>> &coord_to_window_index_map,std::vector<std::pair<Coordinate,int>> &path,ssize_t end_target_position);
float DTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length, std::vector<std::vector<Coordinate>> &window,std::vector<std::pair<Coordinate,int>> &path,ssize_t &end_target_position);
float _fastDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length, int radius,std::vector<std::pair<Coordinate,int>> &path,std::vector<std::vector<Coordinate>> &window, ssize_t &end_target_position);
std::string print_alignment(std::vector<std::pair<Coordinate,int>> &path);
float fastDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length,int radius,ssize_t &end_target_position,std::string &cigar);
} // namespace sigmap

#endif // FAST_DTW_H_
