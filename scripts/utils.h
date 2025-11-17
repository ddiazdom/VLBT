//
// Created by Diaz, Diego on 1.6.2025.
//

#ifndef SR_INDEX_PARSE_PATTERN_H
#define SR_INDEX_PARSE_PATTERN_H

#define MEASURE(query, time_answer, query_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query_answer = query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

#define MEASURE_VOID(query, time_answer, time_unit) \
{\
auto t1 = std::chrono::high_resolution_clock::now();\
query;\
auto t2 = std::chrono::high_resolution_clock::now();\
time_answer += std::chrono::duration_cast<time_unit>( t2 - t1 ).count();\
}

using ulint = uint64_t;
//parse pizza&chilli patterns header:
void pat_error(){
    std::cout << "Error: malformed header in patterns file" << std::endl;
    std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
    exit(0);
}

inline ulint get_n_patterns(std::string header){

    ulint start_pos = header.find("number=");
    if (start_pos == std::string::npos or start_pos+7>=header.size())
        pat_error();

    start_pos += 7;

    ulint end_pos = header.substr(start_pos).find(" ");
    if (end_pos == std::string::npos)
        pat_error();

    ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());
    return n;
}

inline ulint get_patterns_len(std::string header){

    ulint start_pos = header.find("length=");
    if (start_pos == std::string::npos or start_pos+7>=header.size())
        pat_error();

    start_pos += 7;

    ulint end_pos = header.substr(start_pos).find(" ");
    if (end_pos == std::string::npos)
        pat_error();

    ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

    return n;
}

inline std::vector<std::string> file2pat_list(const std::string& pat_file, ulint &n_pats, ulint& pat_len){
    std::ifstream ifs(pat_file);
    std::string header;
    std::getline(ifs, header);
    n_pats = get_n_patterns(header);
    pat_len = get_patterns_len(header);
    std::cout<<"Searching for "<<n_pats<<" patterns of length "<<pat_len<<" each "<<std::endl;
    std::vector<std::string> pat_list(n_pats);
    for(ulint i=0;i<n_pats;++i){
        pat_list[i].reserve(pat_len);
        for(ulint j=0;j<pat_len;++j){
            char c;
            ifs.get(c);
            pat_list[i].push_back(c);
        }
    }
    return pat_list;
}
#endif //SR_INDEX_PARSE_PATTERN_H
