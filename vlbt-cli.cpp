//
// Created by Diaz, Diego on 17.11.2025.
//
#include "CLI11.hpp"
#include "scripts/utils.h"
#include "include/vlbt_build_bwt.h"
#include "include/vlbt_build_sr_index.h"
#include "include/vlbt_sr_index.h"
#include <filesystem>

struct arguments{
    std::string input_file;
    std::string output_file;
    std::string tmp_dir;
    std::string pat_file;
    size_t b_size=4096;
    size_t samp=4;
    VLBT_TYPE dt{};
};

class MyFormatter final : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

template<class dt_type>
void count_int3(const std::string& input_index, const std::string& pat_file, std::string index_name){
    dt_type dt;
    load_from_file(input_index, dt);

    const double bps = double(std::filesystem::file_size(input_index)*8)/double(dt.size());

    uint64_t n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    double acc_time=0;
    size_t acc_count=0;
    std::pair<uint64_t, uint64_t> ans;
    for(auto const& p : pat_list) {
        MEASURE(dt.count(p), acc_time, ans, std::chrono::nanoseconds)
        acc_count+=ans.second-ans.first+1;
    }
    //std::cout<<std::fixed<<std::setprecision(3);
    const double ns_per_pat = acc_time/double(n_pats);
    const double ns_per_occ = acc_time/double(acc_count);
    std::cout<<acc_time<<std::endl;

    std::cout<<"#index_name\tbits_per_sym\tn_pats\tpat_len\tn_occ\tnanosecs/pat\tnanosecs/occ"<<std::endl;
    std::cout<<index_name<<"\t"<<bps<<"\t"<<n_pats<<"\t"<<pat_len<<"\t"<<acc_count<<"\t"<<ns_per_pat<<"\t"<<ns_per_occ<<std::endl;
}

template<size_t b_size>
void count_int2(temp_param_t& tp, const std::string& input_index, const std::string& pat_file) {
    switch (tp.tag) {
        case RLBWT:
            count_int3<vlbt_rlbwt<b_size>>(input_index, pat_file, "RLBWT_"+std::to_string(b_size));
            break;
        case RLBWT_WITH_TOEHOLDS:
            count_int3<vlbt_rlbwt_th<b_size>>(input_index, pat_file, "RLBWT_WITH_TOEHOLDS_"+std::to_string(b_size));
            break;
        case SRI_VALID_AREA:
            count_int3<vlbt_sri_va<b_size, b_size>>(input_index, pat_file, "SRI_VALID_AREA_"+std::to_string(b_size));
            break;
        default:
            exit(1);
    }
}

void count_int(const std::string& input_index, const std::string& pat_file) {

    temp_param_t tp = read_template_param(input_index);
    switch (tp.b_size) {
        case 1024:
            count_int2<1024>(tp, input_index, pat_file);
            break;
        case 4096:
            count_int2<4096>(tp, input_index, pat_file);
            break;
        case 16384:
            count_int2<16384>(tp, input_index, pat_file);
            break;
        case 65536:
            count_int2<65536>(tp, input_index, pat_file);
            break;
        case 262144:
            count_int2<262144>(tp, input_index, pat_file);
            break;
        case 1048576:
            count_int2<1048576>(tp, input_index, pat_file);
            break;
        default:
            exit(1);
    }
}

template<size_t b_size>
void locate_int2(const std::string& input_index, std::string pat_file) {

    vlbt_sri_va<b_size, b_size> sr_index;
    load_from_file(input_index, sr_index);

    uint64_t n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    for(auto const& p : pat_list) {
        std::vector<uint64_t> occ;
        MEASURE(sr_index.locate(p), acc_time, occ, std::chrono::microseconds)
        acc_count+=occ.size();
    }

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"Index type \""<<"SRI_VALID_AREA"<<"\""<<std::endl;
    std::cout<<"\tTotal number of occurrences "<<acc_count<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(n_pats)<<" microsecs/pat"<<std::endl;
    std::cout<<"\t"<<double(acc_time)/double(acc_count)<<" microsecs/occ"<<std::endl;
}

void locate_int(const std::string& input_index, const std::string& pat_file) {
    temp_param_t tp = read_template_param(input_index);
    assert(tp.tag==SRI_VALID_AREA);

    switch (tp.b_size) {
        case 1024:
            locate_int2<1024>(input_index, pat_file);
            break;
        case 4096:
            locate_int2<4096>(input_index, pat_file);
            break;
        case 16384:
            locate_int2<16384>(input_index, pat_file);
            break;
        case 65536:
            locate_int2<65536>(input_index, pat_file);
            break;
        case 262144:
            locate_int2<262144>(input_index, pat_file);
            break;
        case 1048576:
            locate_int2<1048576>(input_index, pat_file);
            break;
        default:
            exit(1);
    }
}

template<uint64_t b_size>
void build_int2(VLBT_TYPE &dt_type, const std::string& input_text, size_t sri_samp_val, const std::filesystem::path& tmp_path, std::string& output_file) {

    std::string bwt_file = input_text+".bwt";
    assert(std::filesystem::exists(bwt_file));

    if (dt_type==RLBWT) {
        std::cout<<"Building the RLBWT with block size "<<b_size<<std::endl;
        vlbt_rlbwt<b_size> bwt;
        build_bwt(bwt, bwt_file, PLAIN, tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("rlbwt_vlt");
        store_to_file(output_file, bwt);
    } else if (dt_type==RLBWT_WITH_TOEHOLDS) {
        std::cout<<"Building the RLBWT with toeholds, using subsampling "<<sri_samp_val<<" and vlbt block size "<<b_size<<std::endl;
        vlbt_rlbwt_th<b_size> bwt_th;
        build_bwt_th<uint64_t>(bwt_th,  bwt_file, PLAIN, sri_samp_val,   tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("rlbwt_th_vlt");
        store_to_file(output_file, bwt_th);
    } else if (dt_type==SRI_VALID_AREA) {
        std::cout<<"Building the sr-index with valid area, using subsampling "<<sri_samp_val<<" and vlbt block size "<<sri_samp_val<<std::endl;
        vlbt_sri_va<b_size, b_size> sri_va;
        build_sr_index<uint64_t>(sri_va, input_text, PLAIN, sri_samp_val, tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("sri_vlt");
        store_to_file(output_file, sri_va);
    } else {
        exit(1);
    }
}

void build_int(VLBT_TYPE& dt_type, const std::string& input_text, const uint64_t b_size,
               const size_t sri_samp_val, const std::filesystem::path &tmp_path, std::string& output_file){
    switch (b_size) {
        case 1024:
            build_int2<1024>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        case 4096:
            build_int2<4096>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        case 16384:
            build_int2<16384>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        case 65536:
            build_int2<65536>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        case 262144:
            build_int2<262144>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        case 1048576:
            build_int2<1048576>(dt_type, input_text, sri_samp_val, tmp_path, output_file);
            break;
        default:
            exit(1);
    }
}

static void parse_app(CLI::App& app, arguments& args){

	const auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(23);
    app.formatter(fmt);

    // Allowed values
    std::map<std::string, int> valid_values{
            {"1024",     1024},
            {"4096",     4096},
            {"16384",   16384},
            {"65536",   65536},
            {"262144", 262144},
            {"1048576",1048576}
    };

    auto * build = app.add_subcommand("build");
    build->add_option("TEXT", args.input_file, "Input file to be indexed")->required();
    build->add_option("-s,--samp", args.samp, "Subsampling parameter (def 4)")->default_val(4);
    build->add_option("-b,--block-size", args.b_size, "Block size (4096)")->required()->transform(CLI::CheckedTransformer(valid_values));
    build->add_option("-d,--dt-type", args.dt, "Data structure to be constructed (0=RLBWT, 1=RLBWT_THLDS, 2=SRI_VAL_AREA)")->required()->check(CLI::Range(0,2));
    build->add_option("-o,--output", args.output_file, "Output file where the index will be stored");
    build->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /os_tmp/vlbt_xxxx)")-> check(CLI::ExistingDirectory);

    auto * count = app.add_subcommand("count");
    count->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    count->add_option("PAT_FILE", args.pat_file, "List of patterns in Pizza&Chilli format")->check(CLI::ExistingFile)->required();

    auto * locate = app.add_subcommand("locate");
    locate->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    locate->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();

    //auto * bkdown = app.add_subcommand("breakdown");
    //bkdown->add_option("INDEX", args.input_file, "Index to be read")->check(CLI::ExistingFile)->required();
    //bkdown->add_option("-i,--index-type", args.index_type, "Subsample r-index variant (0=standard, 1=valid_marks, 2=valid_area)")->required();
    app.require_subcommand(1,1);
}

/*template<class index_type>
void breakdown_int(std::string input_index){
    index_type index;
    sdsl::load_from_file(index, input_index);
    std::vector<std::pair<std::string, size_t>> parts = index.breakdown();
    std::cout<<"Index file: "<<input_index<<std::endl;
    std::cout<<"Subsampling parameter: "<<index.SubsampleRate()<<std::endl;
    size_t acc=0;
    for(auto const& part : parts){
        acc+=part.second;
    }
    for(auto const& part : parts){
        std::cout<<"\t"<<part.first<<": "<<part.second<<" bytes ("<<100*(double)part.second/(double)acc<<"%)"<<std::endl;
    }
    std::cout<<"Total: "<<acc<<" bytes"<<std::endl;
}*/

int main(int argc, char** argv) {

    arguments args;
    CLI::App app("VLBT data structures");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("build")) {
        if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
        build_int(args.dt, args.input_file, args.b_size, args.samp, args.tmp_dir, args.output_file);
    } else if(app.got_subcommand("count")){
        count_int(args.input_file, args.pat_file);
    } else if(app.got_subcommand("locate")){
        locate_int(args.input_file, args.pat_file);
    } /*else if(app.got_subcommand("breakdown")){
        switch (args.index_type) {
            case SRI_INDEX:
                std::cout<<"Index type: sri"<<std::endl;
                breakdown_int<sri::SrIndex<>>(args.input_file);
                break;
            case SRI_VALID_MARKS:
                std::cout<<"Index type: sri_valid_marks"<<std::endl;
                breakdown_int<sri::SrIndexValidMark<>>(args.input_file);
                break;
            case SRI_VALID_AREA:
                std::cout<<"Index type: sri_valid_area"<<std::endl;
                breakdown_int<sri::SrIndexValidArea<>>(args.input_file);
                break;
            default:
                std::cerr<<"Unknown subsample r-index type"<<std::endl;
                exit(1);
        }
    }*/
    else {
        std::cerr<<" Unknown command "<<std::endl;
        exit(1);
    }
    return 0;
}