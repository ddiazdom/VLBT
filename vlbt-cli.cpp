/*
* VLBT – Variable-Length Blocking Trees
 *
 * Copyright (c) 2026 University of Helsinki
 *
 * This file is part of the VLBT software and is distributed under the
 * BSD 3-Clause License. See the LICENSE file for details.
 */
#include "CLI11.hpp"
#include "scripts/utils.h"
#include "include/vlbt_build_bwt.h"
#include "include/vlbt_build_sr_index.h"
#include "include/vlbt_sr_index.h"
#include "include/logger.h"
#include <filesystem>
#include <version.h>//do not delete it (built dynamically to print the program version)

std::string version_string() {
    std::ostringstream out;
    out << PROJECT_NAME << " " << PROJECT_VERSION << "\n";
    out << "commit: " << GIT_COMMIT << "\n";
    out << "built: " << BUILD_DATE << " " << BUILD_TIME << "\n";
    out << "build: " << BUILD_TYPE << "\n";
    return out.str();
}

struct arguments{
    std::string bwt_file;//bwt
    std::string sa_heads_file;//SA heads
    std::string sa_tails_file;//SA tails
    std::string output_file;
    std::string tmp_dir;
    std::string pat_file;
    std::string index_file;
    size_t b_size=4096;
    size_t samp=4;
    VLBT_TYPE dt{};
    log_level log_lvl=log_level::INFO;
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
    const std::string file = std::filesystem::path(input_index).filename();

    if constexpr (dt_type::tag==RLBWT_WITH_TOEHOLDS || dt_type::tag==SRI_VALID_AREA) {
        index_name= index_name+"_s_"+std::to_string(dt.subsampling_value());
    }

    uint64_t n_pats, pat_len;
    const std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    std::pair<uint64_t, uint64_t> ans;
    for(auto const& p : pat_list) {
        MEASURE(dt.count(p), acc_time, ans, std::chrono::nanoseconds)
        acc_count+=ans.second-ans.first+1;
    }
    const double ns_per_pat = double(acc_time)/double(n_pats);
    const double ns_per_occ = double(acc_time)/double(acc_count);

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"#file\tindex_type\tbits_per_sym\tn_pats\tpat_len\tn_occ\tnanosecs/pat\tnanosecs/occ"<<std::endl;
    std::cout<<file<<"\t"<<index_name<<"\t"<<bps<<"\t"<<n_pats<<"\t"<<pat_len<<"\t"<<acc_count<<"\t"<<ns_per_pat<<"\t"<<ns_per_occ<<std::endl;
}

template<size_t b_size>
void count_int2(temp_param_t& tp, const std::string& input_index, const std::string& pat_file) {
    switch (tp.tag) {
        case RLBWT:
            count_int3<vlbt_rlbwt<b_size>>(input_index, pat_file, "RLBWT_b_"+std::to_string(b_size));
            break;
        case RLBWT_WITH_TOEHOLDS:
            count_int3<vlbt_rlbwt_th<b_size>>(input_index, pat_file, "RLBWT_WITH_TOEHOLDS_b_"+std::to_string(b_size));
            break;
        case SRI_VALID_AREA:
            count_int3<vlbt_sri_va<b_size, b_size>>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(b_size));
            break;
        default:
            LOG_ERROR("Unrecognized index type "+std::to_string(tp.tag)+" in "+input_index);
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
            LOG_ERROR("Unsupported block size "+std::to_string(tp.b_size));
            exit(1);
    }
}

template<size_t b_size>
void locate_int2(const std::string& input_index, std::string pat_file, std::string index_name) {

    vlbt_sri_va<b_size, b_size> sr_index;
    load_from_file(input_index, sr_index);

    const double bps = double(std::filesystem::file_size(input_index)*8)/double(sr_index.size());
    const std::string file = std::filesystem::path(input_index).filename();
    index_name= index_name+"_s_"+std::to_string(sr_index.subsampling_value());

    uint64_t n_pats, pat_len;
    std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    for(auto const& p : pat_list) {
        std::vector<uint64_t> occ;
        MEASURE(sr_index.locate(p), acc_time, occ, std::chrono::nanoseconds)
        acc_count+=occ.size();
    }

    const double ns_per_pat = double(acc_time)/double(n_pats);
    const double ns_per_occ = double(acc_time)/double(acc_count);

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"#file\tindex_type\tbits_per_sym\tn_pats\tpat_len\tn_occ\tnanosecs/pat\tnanosecs/occ"<<std::endl;
    std::cout<<file<<"\t"<<index_name<<"\t"<<bps<<"\t"<<n_pats<<"\t"<<pat_len<<"\t"<<acc_count<<"\t"<<ns_per_pat<<"\t"<<ns_per_occ<<std::endl;
}

void locate_int(const std::string& input_index, const std::string& pat_file) {
    temp_param_t tp = read_template_param(input_index);

    if (tp.tag!=SRI_VALID_AREA) {
        LOG_ERROR("Locate requires a SRI_VALID_AREA index, but "+input_index
                  +" has type "+std::to_string(tp.tag));
        exit(1);
    }

    switch (tp.b_size) {
        case 1024:
            locate_int2<1024>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        case 4096:
            locate_int2<4096>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        case 16384:
            locate_int2<16384>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        case 65536:
            locate_int2<65536>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        case 262144:
            locate_int2<262144>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        case 1048576:
            locate_int2<1048576>(input_index, pat_file, "SRI_VALID_AREA_b_"+std::to_string(tp.b_size));
            break;
        default:
            LOG_ERROR("Unsupported block size "+std::to_string(tp.b_size));
            exit(1);
    }
}

template<uint64_t b_size>
void build_int2(VLBT_TYPE &dt_type, const std::string& bwt_file,
                const std::string& sa_heads_file, const std::string& sa_tails_file,
                size_t sri_samp_val, const std::filesystem::path& tmp_path,
                std::string& output_file) {

    if (dt_type==RLBWT) {

        LOG_INFO("Building the RLBWT with block size "+std::to_string(b_size));

        vlbt_rlbwt<b_size> bwt;
        build_bwt(bwt, bwt_file, PLAIN, tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("rlbwt_vlt");
        store_to_file(output_file, bwt);
    } else if (dt_type==RLBWT_WITH_TOEHOLDS) {

        LOG_INFO("Building the RLBWT with toeholds, using subsampling "
                 +std::to_string(sri_samp_val)+" and vlbt block size "
                 +std::to_string(b_size));

        vlbt_rlbwt_th<b_size> bwt_th;
        build_bwt_th<uint64_t>(bwt_th,  bwt_file, sa_heads_file, sa_tails_file, PLAIN, sri_samp_val,   tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("rlbwt_th_vlt");
        store_to_file(output_file, bwt_th);
    } else if (dt_type==SRI_VALID_AREA) {

        LOG_INFO("Building the sr-index with valid area, using subsampling "
                 +std::to_string(sri_samp_val)+" and vlbt block size "+
                 std::to_string(b_size));

        vlbt_sri_va<b_size, b_size> sri_va;
        build_sr_index<uint64_t>(sri_va, bwt_file, sa_heads_file, sa_tails_file, PLAIN, sri_samp_val, tmp_path);
        output_file = std::filesystem::path(output_file).replace_extension("sri_vlt");
        store_to_file(output_file, sri_va);
    } else {
        exit(1);
    }
}

void build_int(VLBT_TYPE& dt_type, const std::string& bwt_file,
               const std::string& sa_heads_file, const std::string& sa_tails_file,
               const uint64_t b_size, const size_t sri_samp_val,
               const std::filesystem::path &tmp_path, std::string& output_file){
    switch (b_size) {
        case 1024:
            build_int2<1024>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        case 4096:
            build_int2<4096>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        case 16384:
            build_int2<16384>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        case 65536:
            build_int2<65536>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        case 262144:
            build_int2<262144>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        case 1048576:
            build_int2<1048576>(dt_type, bwt_file, sa_heads_file, sa_tails_file, sri_samp_val, tmp_path, output_file);
            break;
        default:
            LOG_ERROR("Unsupported block size "+std::to_string(b_size));
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
    build->add_option("BWT", args.bwt_file, "Input BWT to be encoded")->required();
    build->add_option("SAH", args.sa_heads_file, "Run heads of the suffix array");
    build->add_option("SAT", args.sa_tails_file, "Run tails of the suffix array");
    build->add_option("-s,--samp", args.samp, "Subsampling parameter (def. 4)")->default_val(4);
    build->add_option("-b,--block-size", args.b_size, "Block size (def. 4096)")->default_val(4096)->transform(CLI::CheckedTransformer(valid_values));
    build->add_option("-d,--dt-type", args.dt, "Data structure to be constructed (0=RLBWT, 1=RLBWT_THLDS, 2=SRI_VAL_AREA)")->required()->check(CLI::Range(0,2));
    build->add_option("-o,--output", args.output_file, "Output file where the index will be stored");
    build->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /os_tmp/vlbt_xxxx)")->check(CLI::ExistingDirectory)->default_val(std::filesystem::temp_directory_path().string());
    build->add_option("-l,--log-level", args.log_lvl, "Verbosity level (ERROR=0, WARN=1, INFO=2, DEBUG=3, TRACE=4, def. 2)")->check(CLI::Range(0, 4));

    auto * count = app.add_subcommand("count");
    count->add_option("INDEX", args.index_file, "Index file")->check(CLI::ExistingFile)->required();
    count->add_option("PAT_FILE", args.pat_file, "List of patterns in Pizza&Chilli format")->check(CLI::ExistingFile)->required();

    auto * locate = app.add_subcommand("locate");
    locate->add_option("INDEX", args.index_file, "Index file")->check(CLI::ExistingFile)->required();
    locate->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();

    app.require_subcommand(1,1);
    app.set_version_flag("-v,--version", version_string(), "Print the software version and exit");
}

int main(int argc, char** argv) {

    arguments args;
    CLI::App app("VLBT data structures");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    Logger::level = args.log_lvl;

    if(app.got_subcommand("build")) {
        if(args.bwt_file=="-") {//from the stdin
            if (args.dt!=RLBWT) {
                LOG_ERROR("RLBWT is the only data structure that can be built from stdin");
                exit(1);
            }
            if (args.output_file.empty()) {
                LOG_ERROR("An output file must be specified when reading from stdin");
                exit(1);
            }
        } else if(args.output_file.empty()) {
            args.output_file = std::filesystem::path(args.bwt_file).filename();
        }
        build_int(args.dt, args.bwt_file, args.sa_heads_file, args.sa_tails_file, args.b_size, args.samp, args.tmp_dir, args.output_file);
    } else if(app.got_subcommand("count")){
        count_int(args.index_file, args.pat_file);
    } else if(app.got_subcommand("locate")){
        locate_int(args.index_file, args.pat_file);
    } else {
        LOG_ERROR("Unknown command");
        exit(1);
    }
    return 0;
}