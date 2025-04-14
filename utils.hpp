//
// Created by Diaz, Diego on 3.3.2022.
//

#ifndef UTILS_VLBWT_H
#define UTILS_VLBWT_H

#include <iostream>
#include <fstream>

#define INT_CEIL(a,b) (a>0? 1+(a-1)/b : 0)

uint8_t sym_width(unsigned long val){
    if(val==0) return 0;
    return (sizeof(unsigned long)*8) - __builtin_clzl(val);
}

size_t next_power_of_two(unsigned long val){
    uint8_t width = sym_width(val);
    return 1UL<<width;
}

size_t prev_power_of_two(unsigned long val){
    uint8_t width = sym_width(val);
    return 1UL<<(width-1);
}


bool is_power_of_two(unsigned long val){
    return !(val & (val-1));
}

bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s = std::filesystem::file_status{}){
    if(std::filesystem::status_known(s) ? std::filesystem::exists(s) : std::filesystem::exists(p)){
        return true;
    }else{
        return false;
    }
}

std::string random_string(size_t length){
    const std::string characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<size_t> distribution(0, characters.size() - 1);
    std::string random_string;
    random_string.resize(length);
    for(std::size_t i = 0; i < length; i++) {
        random_string[i] = characters[distribution(generator)];
    }
    return random_string;
}

struct tmp_workspace{

    std::string tmp_folder;
    std::string ext;
    bool remove_all;

    explicit tmp_workspace(std::string const& base_folder=std::filesystem::temp_directory_path(),
                           bool rem_all=true,
                           std::string const& prefix="tmp") : remove_all(rem_all) {

        std::string tmp_path = std::filesystem::canonical(std::filesystem::path(base_folder)) / std::string(prefix+".XXXXXX");
        char temp[200] = {0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size() + 1] = '\0';
        auto res = mkdtemp(temp);
        if (res == nullptr) {
            std::cout << "Error trying to create a temporal folder" << std::endl;
        }
        tmp_folder = std::string(temp);
        ext = random_string(3);
    }

    explicit tmp_workspace(std::string folder,
                           std::string ext_,
                           bool rem_all) : tmp_folder(folder),
                                           ext(ext_),
                                           remove_all(rem_all){
        assert(file_exists(tmp_folder));
    }

    [[nodiscard]] std::string get_file(std::string const& prefix) const {
        return std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
    }

    void remove_file(std::string const& prefix) const {
        std::filesystem::path file =  std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
        bool res = remove(file);
        if(!res){
            std::cout<<"Error trying to remove "<<file<<std::endl;
            exit(1);
        }
    }

    ~tmp_workspace(){
        if(remove_all){
            std::filesystem::remove_all(tmp_folder);
        }
    }

    [[nodiscard]] std::string folder() const{
        return tmp_folder;
    }
};

size_t round_to_power_of_two(unsigned long val){
    if(is_power_of_two(val)){
        return val;
    }else{
        return next_power_of_two(val);
    }
}

template<class vector_t>
size_t basic_store_vector_to_file(std::string const& file, vector_t& vector){
    std::ofstream ofs(file, std::ios::binary);
    ofs.write((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*vector.size()));
    ofs.close();
    return sizeof(typename vector_t::value_type)*vector.size();
}

template<class vector_t>
size_t basic_load_vector_from_file(std::string const& file, vector_t& vector){
    std::ifstream ifs(file, std::ios::binary);
    ifs.seekg(0, std::ios::end);
    std::streamsize n_bytes = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    size_t type_bytes = sizeof(typename vector_t::value_type);
    vector.resize(n_bytes/type_bytes);
    ifs.read((char *)vector.data(), n_bytes);
    ifs.close();
    return n_bytes;
}

template<class data_type>
void load_from_file(std::string const& file, data_type& dt){
    std::ifstream ifs(file, std::ios::binary);
    dt.load(ifs);
    ifs.close();
}

template<class data_type>
size_t store_to_file(std::string const& file, data_type& dt){
    std::ofstream ofs(file, std::ios::binary);
    size_t written_bytes = dt.serialize(ofs);
    ofs.close();
    return written_bytes;
}

template<uint8_t width>
static inline bool aligned(size_t bit_pos){
    return (bit_pos % width)==0;
}

template<class vector_t>
size_t serialize_plain_vector(std::ostream& ofs, vector_t& vector){
    size_t n = vector.size();
    ofs.write((char *)&n, sizeof(n));
    ofs.write((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
    return sizeof(n)+ sizeof(typename vector_t::value_type)*n;
}

template<class size_type>
size_t serialize_raw_vector(std::ostream& ofs, size_type * vector, size_t len){
    ofs.write((char *)&len, sizeof(len));
    ofs.write((char *)vector, (std::streamsize)(sizeof(size_type)*len));
    return sizeof(len)+ sizeof(size_type)*len;
}

template<class val_type>
size_t serialize_elm(std::ostream& ofs, val_type value){
    ofs.write((char *)&value, sizeof(val_type));
    return sizeof(val_type);
}

template<class vector_t>
void load_plain_vector(std::istream& ifs, vector_t& vector){
    size_t n=0;
    ifs.read((char *)&n, sizeof(n));
    vector.resize(n);
    ifs.read((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
}

template<class size_type, class len_type>
void load_raw_vector(std::istream& ifs, size_type*& vector, len_type& len){
    ifs.read((char *)&len, sizeof(len_type));
    if(vector== nullptr){
        vector = (size_type *) malloc(sizeof(size_type)*len);
    }else{
        vector = (size_type *) realloc(vector, sizeof(size_type)*len);
    }
    ifs.read((char *)vector, (std::streamsize)(sizeof(size_type)*len));
}

template<class val_type>
void load_elm(std::istream& ifs, val_type& value){
    ifs.read((char *)&value, sizeof(val_type));
}

template<class vector_type>
void store_pl_vector(std::string const& file, vector_type& vector){
    std::ofstream ofs(file, std::ios::binary);
    serialize_plain_vector<vector_type>(ofs, vector);
    ofs.close();
}

template<class vector_type>
void load_pl_vector(std::string const& file, vector_type& vector){
    std::ifstream ifs(file, std::ios::binary);
    load_plain_vector<vector_type>(ifs, vector);
    ifs.close();
}
#endif //UTILS_VLBWT_H
