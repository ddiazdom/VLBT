#include "CLI11.hpp"
#include "include/sr-index/sr_index.h"
#include "include/sr-index/construct.h"
#include "include/sr-index/config.h"
#include "sri_cli_utils.h"

#include <filesystem>
#include <random>

// Helper: generate a random hex string for unique names
std::string random_hex(std::size_t length = 16) {
    static std::mt19937_64 rng{std::random_device{}()};
    static std::uniform_int_distribution<uint64_t> dist;

    std::ostringstream oss;
    while (oss.tellp() < static_cast<std::streamoff>(length)) {
        oss << std::hex << std::setw(16) << std::setfill('0') << dist(rng);
    }
    return oss.str().substr(0, length);
}

namespace fs = std::filesystem;
std::string create_tmp_dir(std::string prefix=""){
    fs::path temp_root;
    if (prefix.empty()) {
        temp_root = fs::temp_directory_path();
    }else {
        temp_root = prefix;
    }

    fs::path temp_dir = temp_root / fs::path("sri_"+random_hex(12));
    while (fs::exists(temp_dir)) {
        temp_dir = temp_root / fs::path("sri_"+random_hex(12));
    }
    fs::create_directories(temp_dir);

    return std::filesystem::absolute(temp_dir).lexically_normal().string();
}

enum SRI_TYPE{
    SRI_INDEX=0,
    SRI_VALID_MARKS=1,
    SRI_VALID_AREA=2,
};

struct arguments{
    std::string input_file;
    std::string output_file;
    std::string tmp_dir="";
    std::string pat_file;
    size_t n_threads{};
    size_t ssamp=4;
    sri::SAAlgo sa_algo = sri::SDSL_LIBDIVSUFSORT;
    SRI_TYPE index_type = SRI_VALID_AREA;
    std::string bigbwt_pref;
    size_t bytes_sa=5;
};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

template<class index_type>
void test_count(std::string input_file, std::string& pat_file, std::string index_name){

    index_type index;
    sdsl::load_from_file(index, input_file);

    const double bps = double(std::filesystem::file_size(input_file)*8)/double(index.sizeSequence());
    const std::string file = std::filesystem::path(input_file).filename();
    index_name=index_name+"_s_"+std::to_string(index.SubsampleRate());

    uint64_t n_pats, pat_len;
    const std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    std::pair<size_t, size_t> ans;
    for(auto const& p : pat_list) {
        MEASURE(index.Count(p), acc_time, ans, std::chrono::nanoseconds)
        acc_count+=ans.second-ans.first+1;
    }

    //std::cout<<std::fixed<<std::setprecision(3);
    //std::cout<<"Index type \""<<index_name<<"\""<<std::endl;
    //std::cout<<"\tTotal number of occurrences "<<acc_count<<std::endl;
    //std::cout<<"\t"<<double(acc_time)/double(n_pats)<<" nanosecs/pat"<<std::endl;
    //std::cout<<"\t"<<double(acc_time)/double(acc_count)<<" nanosecs/occ"<<std::endl;

    const double ns_per_pat = double(acc_time)/double(n_pats);
    const double ns_per_occ = double(acc_time)/double(acc_count);

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"#file\tindex_type\tbits_per_sym\tn_pats\tpat_len\tn_occ\tnanosecs/pat\tnanosecs/occ"<<std::endl;
    std::cout<<file<<"\t"<<index_name<<"\t"<<bps<<"\t"<<n_pats<<"\t"<<pat_len<<"\t"<<acc_count<<"\t"<<ns_per_pat<<"\t"<<ns_per_occ<<std::endl;
}

template<class index_type>
void test_locate(std::string input_file, std::string& pat_file, std::string index_name){

    index_type index;
    sdsl::load_from_file(index, input_file);

    const double bps = double(std::filesystem::file_size(input_file)*8)/double(index.sizeSequence());
    const std::string file = std::filesystem::path(input_file).filename();
    index_name=index_name+"_s_"+std::to_string(index.SubsampleRate());

    uint64_t n_pats, pat_len;
    const std::vector<std::string> pat_list = file2pat_list(pat_file, n_pats, pat_len);

    size_t acc_time=0;
    size_t acc_count=0;
    for(auto const& p : pat_list) {
        std::vector<size_t> occ;
        MEASURE(index.Locate(p), acc_time, occ, std::chrono::nanoseconds)
        acc_count+=occ.size();
    }

    const double ns_per_pat = double(acc_time)/double(n_pats);
    const double ns_per_occ = double(acc_time)/double(acc_count);

    std::cout<<std::fixed<<std::setprecision(3);
    std::cout<<"#file\tindex_type\tbits_per_sym\tn_pats\tpat_len\tn_occ\tnanosecs/pat\tnanosecs/occ"<<std::endl;
    std::cout<<file<<"\t"<<index_name<<"\t"<<bps<<"\t"<<n_pats<<"\t"<<pat_len<<"\t"<<acc_count<<"\t"<<ns_per_pat<<"\t"<<ns_per_occ<<std::endl;

    //std::cout<<std::fixed<<std::setprecision(3);
    //std::cout<<"Index type \""<<index_name<<"\""<<std::endl;
    //std::cout<<"\tTotal number of occurrences "<<acc_count<<std::endl;
    //std::cout<<"\t"<<double(acc_time)/double(n_pats)<<" microsecs/pat"<<std::endl;
    //std::cout<<"\t"<<double(acc_time)/double(acc_count)<<" microsecs/occ"<<std::endl;
}

static void parse_app(CLI::App& app, struct arguments& args){
    
	auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(23);
    app.formatter(fmt);

    auto * build = app.add_subcommand("build");
    build->add_option("-s,--ssamp", args.ssamp, "Subsampling parameter (def 4)")->default_val(4);
    build->add_option("-i,--index-type", args.index_type, "Subsample r-index variant to be constructed (0=standard, 1=valid_marks, 2=valid_area [def=2])")->default_val(SRI_VALID_AREA)->check(CLI::Range(0,2));
    build->add_option("-t,--threads", args.n_threads, "Maximum number of working threads")->default_val(1);
    build->add_option("-o,--output", args.output_file, "Output file where the index will be stored");
    build->add_option("-T,--tmp", args.tmp_dir, "Temporary folder (def. /os_tmp/sri_xxxx)")-> check(CLI::ExistingDirectory);
    auto *build_algo = build->add_option("-a,--sa-algorithm", args.sa_algo, "Algorithm for computing the SA (0=LIBDIVSUFSORT, 1=SE_SAIS, 2=BIG_BWT [def=0])")->default_val(LIBDIVSUFSORT)->check(CLI::Range(0,2));

    auto * group_option = build->add_option_group("Source of the index components (one of the two is mandatory):");
    auto *text = group_option->add_option("-t,--text", args.input_file, "Input TEXT to be indexed")->check(CLI::ExistingFile);

    auto *bwt_pref = group_option->add_option("-b,--bigbwt-pref", args.bigbwt_pref, "BigBWT prefix containing the already computed BWT and SA samples")->expected(1)->check(CLI::ExistingFile);

    group_option->require_option(1,2);
    bwt_pref->excludes(text);
    bwt_pref->excludes(build_algo);

    build_algo->needs(text);

    auto * count = app.add_subcommand("count");
    count->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    count->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();
    count->add_option("-i,--index-type", args.index_type, "Subsample r-index variant (0=standard, 1=valid_marks, 2=valid_area)")->required();

    auto * locate = app.add_subcommand("locate");
    locate->add_option("INDEX", args.input_file, "Index file")->check(CLI::ExistingFile)->required();
    locate->add_option("PAT_FILE", args.pat_file, "List of patterns")->check(CLI::ExistingFile)->required();
    locate->add_option("-i,--index-type", args.index_type, "Subsample r-index variant (0=standard, 1=valid_marks, 2=valid_area)")->required();

    auto * bkdown = app.add_subcommand("breakdown");
    bkdown->add_option("INDEX", args.input_file, "Index to be read")->check(CLI::ExistingFile)->required();
    bkdown->add_option("-i,--index-type", args.index_type, "Subsample r-index variant (0=standard, 1=valid_marks, 2=valid_area)")->required();

    app.require_subcommand(1,1);
}

template<class index_type>
void build_int(std::string input_text, size_t ssamp_val, std::filesystem::path tmp_path, sri::SAAlgo sa_algo, std::string& output_file){
    index_type index(ssamp_val);
    sri::Config config(input_text, tmp_path, sa_algo);
    sri::construct(index, input_text, config);
    sdsl::store_to_file(index, output_file);
}

template<class index_type>
void build_from_bigbwt(std::string bigbwt_pref, size_t ssamp_val, std::filesystem::path tmp_path, std::string& output_file){
    index_type index(ssamp_val);
    sri::Config config(bigbwt_pref, tmp_path, sri::SAAlgo::BIG_BWT);
    sri::construct(index, bigbwt_pref, config);
    sdsl::store_to_file(index, output_file);
}

template<class index_type>
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
}

int main(int argc, char** argv) {

    arguments args;
    CLI::App app("Subsample r-index");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("build")) {

        std::string tmp_dir = create_tmp_dir(args.tmp_dir);
        std::cout<<"Temporary folder: "<<tmp_dir<<std::endl;

        if (!args.input_file.empty()) {
            assert(args.bigbwt_pref.empty());
            if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
            std::cout<<"Building the subsample r-index for "<<args.input_file<<std::endl;
            std::cout<<"Subsampling parameter: "<<args.ssamp<<std::endl;
            switch (args.index_type) {
                case SRI_INDEX:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri");
                    build_int<sri::SrIndex<>>(args.input_file, args.ssamp, tmp_dir, args.sa_algo, args.output_file);
                    break;
                case SRI_VALID_MARKS:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri_vm");
                    build_int<sri::SrIndexValidMark<>>(args.input_file, args.ssamp, tmp_dir, args.sa_algo, args.output_file);
                    break;
                case SRI_VALID_AREA:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri_va");
                    build_int<sri::SrIndexValidArea<>>(args.input_file, args.ssamp, tmp_dir, args.sa_algo, args.output_file);
                    break;
                default:
                    std::cerr<<"Unknown subsample r-index type"<<std::endl;
                    exit(1);
            }
            std::cout<<"The output index was stored in "<<args.output_file<<std::endl;
        } else if (!args.bigbwt_pref.empty()) {
            assert(args.input_file.empty());

            if(args.output_file.empty()) args.output_file = std::filesystem::path(args.bigbwt_pref).filename();

            std::cout<<"Building the subsample r-index from the precomputed BWT/SA elements in "<<args.bigbwt_pref<<std::endl;
            std::cout<<"Subsampling parameter: "<<args.ssamp<<std::endl;
            switch (args.index_type) {
                case SRI_INDEX:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri");
                    build_from_bigbwt<sri::SrIndex<>>(args.bigbwt_pref, args.ssamp, tmp_dir, args.output_file);
                    break;
                case SRI_VALID_MARKS:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri_vm");
                    build_from_bigbwt<sri::SrIndexValidMark<>>(args.bigbwt_pref, args.ssamp, tmp_dir, args.output_file);
                    break;
                case SRI_VALID_AREA:
                    args.output_file = std::filesystem::path(args.output_file).replace_extension("sri_va");
                    build_from_bigbwt<sri::SrIndexValidArea<>>(args.bigbwt_pref, args.ssamp, tmp_dir, args.output_file);
                    break;
                default:
                    std::cerr<<"Unknown subsample r-index type"<<std::endl;
                    exit(1);
            }
            std::uintmax_t removed = fs::remove_all(tmp_dir);
            std::cout<<"The output index was stored in "<<args.output_file<<std::endl;
        }

    } else if(app.got_subcommand("count")){
        switch (args.index_type) {
            case SRI_INDEX:
                test_count<sri::SrIndex<>>(args.input_file, args.pat_file, "sri");
                break;
            case SRI_VALID_MARKS:
                test_count<sri::SrIndexValidMark<>>(args.input_file, args.pat_file, "sri_valid_marks");
                break;
            case SRI_VALID_AREA:
                test_count<sri::SrIndexValidArea<>>(args.input_file, args.pat_file, "sri_valid_area");
                break;
            default:
                std::cerr<<"Unknown subsample r-index type"<<std::endl;
                exit(1);
        }
    } else if(app.got_subcommand("locate")){
        switch (args.index_type) {
            case SRI_INDEX:
                test_locate<sri::SrIndex<>>(args.input_file, args.pat_file, "sri");
                break;
            case SRI_VALID_MARKS:
                test_locate<sri::SrIndexValidMark<>>(args.input_file, args.pat_file, "sri_valid_marks");
                break;
            case SRI_VALID_AREA:
                test_locate<sri::SrIndexValidArea<>>(args.input_file, args.pat_file, "sri_valid_area");
                break;
            default:
                std::cerr<<"Unknown subsample r-index type"<<std::endl;
                exit(1);
        }
    } else if(app.got_subcommand("breakdown")){
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
    } else {
        std::cerr<<" Unknown command "<<std::endl;
        exit(1);
    }
    return 0;
}
