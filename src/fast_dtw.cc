#include "fast_dtw.h"
namespace sigmap {
void reduce_by_half(const float *signal, size_t length, float *signal_reduced, size_t &reduced_length){
    reduced_length = 0;
    for (size_t i = 0; i < length - length % 2; i += 2){
        signal_reduced[reduced_length] = (signal[i] + signal[i + 1]) / 2;
        reduced_length++;
    }
}

void expand_window(std::vector<std::pair<Coordinate,int>> &path,size_t target_length,size_t query_length,int radius,std::vector<std::vector<Coordinate>> &window) {
    std::set<Coordinate> path_set;
    for (size_t i = 0;i<path.size();i++){
        for (int j=-radius;j<=radius;j++){
            for (int k=-radius;k<=radius;k++){
                ssize_t new_target_coord = path[i].first.target_coord+j;
                ssize_t new_query_coord = path[i].first.query_coord+k;
                if ((new_target_coord>=0)&(new_query_coord>=0)&(new_target_coord<target_length)&(new_query_coord<query_length)){
                    path_set.insert(Coordinate(new_target_coord,new_query_coord));
                }
            }
        }
    }
    std::set<Coordinate> window_set;
    for (std::set<Coordinate>::iterator it=path_set.begin() ;it!=path_set.end();it++){
        for (int x=0;x<2;x++){
            for (int y=0;y<2;y++){
                ssize_t new_target_coord = it->target_coord*2+x;
                ssize_t new_query_coord = it->query_coord*2+y;
                if ((new_target_coord>=0)&(new_query_coord>=0)&(new_target_coord<target_length)&(new_query_coord<query_length)){
                    window_set.insert(Coordinate(new_target_coord,new_query_coord));
                }
             }
        }
    }


    //each vector in window represents a base 
    std::vector<std::vector<Coordinate>>().swap(window);
    size_t last_target_coord = window_set.begin()->target_coord;
    window.push_back(std::vector<Coordinate>());
    for (std::set<Coordinate>::iterator it=window_set.begin() ;it!=window_set.end();it++){
        if (last_target_coord!=it->target_coord){
            window.push_back(std::vector<Coordinate>());
            last_target_coord = it->target_coord;
        }
        window.back().push_back(Coordinate(it->target_coord,it->query_coord));
    }
}

void generate_path(std::vector<std::vector<int>> &path_matrix,std::vector<std::vector<Coordinate>> &window,std::map<Coordinate,std::pair<size_t,size_t>> &coord_to_window_index_map,std::vector<std::pair<Coordinate,int>> &path,ssize_t end_target_position){
    //trace back
    size_t row = end_target_position;
    size_t col = window[end_target_position].size()-1;
    //clear the path
    std::vector<std::pair<Coordinate,int>>().swap(path);
    Coordinate coord = window[row][col];
    while (coord.query_coord != 0) {
        coord = window[row][col];
        path.push_back(std::make_pair(Coordinate(coord.target_coord,coord.query_coord),path_matrix[row][col]));
        //coordinate shift when trace back
        int query_shift[] =  {-1,-1,-1,0};
        int target_shift[] = {-1,0,0,-1};
        coord.query_coord += query_shift[path_matrix[row][col]];
        coord.target_coord += target_shift[path_matrix[row][col]];
        std::pair<size_t,size_t> index = coord_to_window_index_map[coord];
        row = index.first;
        col = index.second;
    }
    path.push_back(std::make_pair(window[row][col],path_matrix[row][col])); 
    std::reverse(path.begin(),path.end());
}

float DTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length,  std::vector<std::vector<Coordinate>> &window,std::vector<std::pair<Coordinate,int>> &path,ssize_t &end_target_position){   
    float min_dtw_distance = std::numeric_limits<float>::max();
    float skip_target_cost = 2;
    float skip_query_cost = 2;
    std::map<Coordinate,std::pair<size_t,size_t>> coord_to_window_index_map;
    if ((window.empty())){
        for (size_t i = 0; i < target_length; i++){
            window.push_back(std::vector<Coordinate>());
            for (size_t j = 0; j < query_length; j++){
                window.back().push_back(Coordinate(i,j));
            }
        }
    }

    //path_matrix: 0 -> one-to-one map; 1 -> one base to multi signal; 2 -> skip one signal; 3 -> skip one base;
    //initialize path_matrix and coord_to_window_index_map
    std::vector<std::vector<int>> path_matrix(window.size(),std::vector<int>());
    for (size_t i = 0; i < window.size(); i++){
        path_matrix[i]= std::vector<int>(window[i].size(),0);
        for (size_t j = 0; j < window[i].size(); j++){
            coord_to_window_index_map[window[i][j]] = std::make_pair(i,j);
        }
    }
    std::vector<float> previous_row(query_length + 1, std::numeric_limits<float>::max());
    std::vector<float> current_row(query_length + 1, std::numeric_limits<float>::max());
    previous_row[0] = 0;
    end_target_position = -1;
    size_t target_position;
    size_t query_position;
    for (size_t i = 0; i < window.size(); ++i){
        // iterate to a new row (new signal on reference)
        current_row[0] = 0;
        for (size_t j = 0; j < window[i].size(); ++j){
            target_position = window[i][j].target_coord + 1;
            query_position = window[i][j].query_coord + 1;
            float cost = std::abs(target_signal[target_position - 1] - query_signal[query_position - 1]);
            std::vector<float> candidates;
            std::vector<int> candidates_path_direction;
            //for the start position
            candidates = {previous_row[query_position-1]+cost,current_row[query_position-1]+cost,current_row[query_position-1]+skip_query_cost,previous_row[query_position]+skip_target_cost};
            candidates_path_direction = {0,1,2,3};
            //  if ((target_position==1)&&(query_position==1)){
            //     candidates = {cost,skip_query_cost,skip_target_cost};
            //     candidates_path_direction = {0,2,3};
            // } else if ((j>0) && (path_matrix[i][j-1]==3)){
            //     // if skip this base in the last step, one base to multi signal is not allowed
            //     candidates = {previous_row[query_position-1]+cost,current_row[query_position-1]+skip_query_cost,previous_row[query_position]+skip_target_cost};
            //     candidates_path_direction = {0,2,3};
            // } else {
            //     candidates = {previous_row[query_position-1]+cost,current_row[query_position-1]+cost,current_row[query_position-1]+skip_query_cost,previous_row[query_position]+skip_target_cost};
            //     candidates_path_direction = {0,1,2,3};
            // }
            int argmin_candidates = std::min_element( candidates.begin(),candidates.end()) - candidates.begin();
            current_row[query_position] = candidates[argmin_candidates];
            path_matrix[i][j] = candidates_path_direction[argmin_candidates];
        }
         if ((query_position == query_length) && (current_row[query_position] < min_dtw_distance)) {
            min_dtw_distance = current_row[query_position];
            end_target_position = i;
        }
        current_row.swap(previous_row);
        std::vector<float>(query_length + 1, std::numeric_limits<float>::max()).swap(current_row);
    }
    generate_path(path_matrix,window,coord_to_window_index_map,path,end_target_position);
    end_target_position = window[end_target_position][0].target_coord;
    return min_dtw_distance;
}

float _fastDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length, int radius,std::vector<std::pair<Coordinate,int>> &path,std::vector<std::vector<Coordinate>> &window, ssize_t &end_target_position){
    int min_time_size = radius + 2;
    
     std::vector<std::vector<Coordinate>>().swap(window);
    if ((target_length < min_time_size) || (query_length < min_time_size)){
        return DTW(target_signal, target_length, query_signal, query_length,window,path,end_target_position);
    }
    size_t target_reduced_length = 0;
    size_t query_reduced_length = 0;
    float target_reduced[target_length / 2 + 1];
    float query_reduced[query_length / 2 + 1];
    reduce_by_half(target_signal, target_length, target_reduced, target_reduced_length);
    reduce_by_half(query_signal, query_length, query_reduced, query_reduced_length);
    float dist = _fastDTW(target_reduced,target_reduced_length, query_reduced,query_reduced_length,radius,path,window,end_target_position);
    expand_window(path, target_length,query_length, radius,window);
    return DTW(target_signal, target_length, query_signal, query_length,window,path,end_target_position);
}

std::string print_alignment(std::vector<std::pair<Coordinate,int>> &path){
    //path[i].second: 0 -> one-to-one map; 1 -> one base to multi signal; 2 -> skip one signal; 3 -> skip one base;
    std::string print_flags = "MMID";
    if (path.empty()){
        return "";
    }
    std::vector<std::string> per_base_cigar;
    size_t num_same_continuous_flag = 1;
    int last_flag;
    // for the first base to signal alignment
    if (path[0].second == 3){
        per_base_cigar.push_back("1D");
        last_flag = 3;
    } else {
        last_flag = path[0].second == 0 ? 1 : 2;
    }
    per_base_cigar.push_back("");
    
    for (size_t i = 1; i < path.size(); ++i){
        int flag = path[i].second;
        // if within one base,counting the continous same flag
        if ((flag == 1) || (flag == 2)){
            if (last_flag == flag){
                num_same_continuous_flag ++;
            } else {
                per_base_cigar.back() += std::to_string(num_same_continuous_flag) + print_flags[last_flag];
                num_same_continuous_flag = 1;
                last_flag = flag;
            }     
        } else {
            // if not within one base, finish the last base and adding new base
            per_base_cigar.back() += std::to_string(num_same_continuous_flag) + print_flags[last_flag];
            if (flag == 0){
                last_flag = 1;
            } else if (flag == 3){
                // if skip one base
                last_flag = 3;
            }
            if (i != path.size() - 1){
                per_base_cigar.push_back("");
                num_same_continuous_flag = 1;
            }
        }   
    }
    std::string cigar = "";
    for (size_t i = 0; i < per_base_cigar.size(); ++i){
        cigar += "(" + per_base_cigar[i] + ")";
    }
    return cigar;
}

float fastDTW(const float *target_signal, size_t target_length, const float *query_signal, size_t query_length,int radius,ssize_t& end_target_position,std::string &cigar){
    double real_start_time = GetRealTime();
    std::vector<std::pair<Coordinate,int>> path;
    std::vector<std::vector<Coordinate>> window;
    float distance = _fastDTW(target_signal,target_length,query_signal,query_length,radius,path,window,end_target_position);
    std::cerr << "Finished sDTW in "<< GetRealTime() - real_start_time << ", target length: " << target_length << ", query length: " << query_length << "\n";
    std::cerr<< path[0].first.target_coord<<","<<path[0].first.query_coord<<'\n';
    cigar = print_alignment(path);
    return distance;
}

} // namespace sigmap
