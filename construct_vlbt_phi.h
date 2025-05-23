//
// Created by Diaz, Diego on 17.4.2025.
//

#ifndef VLBT_CONSTRUCT_VLBT_PHI
#define VLBT_CONSTRUCT_VLBT_PHI

#include "vlbt_phi.h"
#ifdef __linux__
#include <malloc.h>
#endif

//statistics about the data structure
template<class phi_dt_type>
struct phi_stat_collector{
    uint64_t rpl_freq[phi_dt_type::max_block_runs+1]={0};//number of runs in a leaf
    uint64_t leaf_depth_freq[20]={0};//the depth of each leaf
    uint64_t leaf_enc_freq[20]={0};//encoding of each leaf
    uint64_t children_freq[100]={0};//children frequency = how many nodes with 1,2,...,x children
    uint64_t header_overhead=0;//number of bits used by the headers of the nodes
    uint64_t runs_overhead=0;
    uint64_t sym_overhead=0;
    uint64_t len_overhead=0;
    uint64_t tree_pointers_overhead=0;
    uint64_t trees_overhead=0;
    uint64_t max_n_blocks=0;
};

template<class phi_dt_type, class size_type>
struct phi_node {//state of the compression

    typedef std::pair<size_type, size_type> run_type;
    typedef std::vector<run_type> block_type;

    size_t bk_len=0;//length of the active block
    size_t bk_id=0;//id of the active block
    size_t acc_runs=0;//sum of the runs in the active blocks
    size_t n_children=0;//number of children of this node
    size_t node_n_bits=0;//number of bits required for the subtree rooted under this node
    size_t consumed_syms=0;//number of symbols under this node that have been scanned so far
    size_t syms_before=0;//number of symbols in the text preceding this node
    size_t child_rank=0;//this node is the child_rank of its parent
    size_t cov_symbols=0;//how many symbols of the input does this node cover?
    size_t child_mark_acc=0;//to mark which children were collapsed
    size_t node_offset=0;//sum of symbols covered by the left siblings of this node
    size_t len_enc=0;//code for the encoding of the run lengths

    bool lm_tree_branch=true;//true if this node in the leftmost branch of its tree
    bool rm_tree_branch=true;//true if this node in the rightmost branch of its tree
    bool lm_child=false;//true if this node is the leftmost child of its parent
    bool rm_child=false;//true if this node is the rightmost child of its parent
    bool leaf=false;//true if the node is a leaf

    const size_t lvl;//level of the subtree
    const size_t b_size;//block size for the level
    const size_t s_factor = phi_dt_type::scale_factor;//shrinking factor for further subdivision
    const size_t b_runs = phi_dt_type::max_block_runs;//maximum number of runs in a sequence of blocks
    const uint8_t run_width = phi_dt_type::run_width;//number of bits we require to encode symbols in the range [0..b_runs]

    phi_node *tmp_node = nullptr;

    stream_type buffer;
    stream_type children_buffer;

    phi_dt_type& phi_rep; //data structure encoding the representation
    std::vector<block_type> active_blocks; //run-length compressed blocks conforming a tree node
    std::vector<bool> child_marks;//int a block of size b_size, it is a bit vector B[1..s_factor] that marks the original blocks of size b_size/s_factor in the collapsed blocks
    std::vector<uint64_t> block_ptr;//pointers (byte offsets) to the node's children

    //this information is for the tree
    std::ostream *ofs= nullptr;
    std::ifstream *ifs= nullptr;
    std::vector<uint64_t> tree_offset;//number of symbols in the text before each tree

    //a struct to collect statistics about the data structure
    phi_stat_collector<phi_dt_type>& stats;

    explicit phi_node(size_t _lvl, size_t _b_size, phi_dt_type& _phi_rep, phi_stat_collector<phi_dt_type>& st):
            lvl(_lvl),
            b_size(_b_size),
            phi_rep(_phi_rep),
            active_blocks(b_runs),
            stats(st){
        if(lvl==0){
            //number of blocks in the tree representation
            size_t n_blocks = INT_CEIL(phi_rep.tot_syms, b_size);
            block_ptr.resize(n_blocks+2);
            block_ptr[0] = 0;
            tree_offset = std::vector<uint64_t>(n_blocks+2, 0);
        } else{
            block_ptr.resize(s_factor+1);
            child_marks.resize(s_factor+1);
        }
    }

    inline void process_block_seq() {
        if(bk_id==1){//only one block in the sequence and it exceeds the limit of runs
            create_node<INTERNAL>(1);//recursive partitioning
            active_blocks[0].clear();
            bk_id=0;
            acc_runs=0;
        } else {//multiple blocks in the block sequence
            create_node<LEAF>(bk_id-1);//we do not consider the last block in the sequence
            active_blocks[0].swap(active_blocks[bk_id-1]);//move the last block to the beginning of the sequence
            for(size_t i=1;i<bk_id;i++){//clear the other blocks in the sequence
                active_blocks[i].clear();
            }

            //new block sequence
            acc_runs = active_blocks[0].size();
            bk_id=1;
            if(active_blocks[0].size()>b_runs){//process the remaining block if it exceeds the limit of runs
                create_node<INTERNAL>(1);
                active_blocks[0].clear();
                bk_id=0;
                acc_runs=0;
            }
        }
    }

    inline void destroy(){
        buffer.destroy();
        children_buffer.destroy();
        for(auto bk : active_blocks){
            destroy_vector(bk);
        }
        destroy_vector(block_ptr);
        destroy_vector(tree_offset);
        destroy_vector(child_marks);
#ifdef __linux__
        malloc_trim(0);
#endif
    }

    inline void finish_run_scan(){
        //handle the last sequence of blocks
        assert(bk_len<=b_size);
        acc_runs+=active_blocks[bk_id].size();
        bk_id+=!active_blocks[bk_id].empty();

        if(acc_runs>=b_runs){//the last block sequence has more than the maximum number of allowed runs
            process_block_seq();
        }

        assert(acc_runs<b_runs);
        if(acc_runs>0){//the block sequence still has some runs left
            create_node<LEAF>(bk_id);
            for(size_t i=0;i<bk_id;i++){//clear the other blocks in the sequence
                active_blocks[i].clear();
            }
            bk_id=0;
            acc_runs=0;
        }
        bk_len=0;
        assert(aligned<8>(node_n_bits));
    }

    inline void process_run(const size_type& sym, size_type& len) {
        while(len>0){
            if((bk_len+len)<b_size){
                //the run fits the block size
                assert(len<=b_size);
                bk_len+=len;
                active_blocks[bk_id].emplace_back(sym, len);
                len=0;
            } else {// we complete a new block
                //last run of the active block
                size_t split_run_len = b_size-bk_len;
                active_blocks[bk_id].emplace_back(sym, split_run_len);
                bk_len+=split_run_len;
                assert(split_run_len>0 && bk_len==b_size);
                acc_runs+=active_blocks[bk_id].size();
                bk_id++;

                //we compete a new block sequence
                if(acc_runs>=b_runs){
                    process_block_seq();
                }
                len -=split_run_len;
                bk_len=0;
            }
        }
    }

    inline void finish_int_node(){

        assert(n_children<=s_factor);
        assert(lvl>0);

        //HEADER DESCRIPTION:
        //1 bit to indicate it is an internal node
        //off_width bits store the offset of this node
        //s bits to indicate which children were collapsed
        //n_children pointers to the children

        //===some preliminary information
        //width in bits for the offset information
        size_t off_width;
        if(lvl==1) {
            //the root of the tree
            off_width = sym_width(phi_rep.tot_syms);
        } else {
            //internal node that is not the root
            //bsize*s_factor*s_factor is the block size of the parent
            off_width = sym_width(b_size*s_factor*s_factor);
            assert(node_offset<(b_size*s_factor*s_factor));
        }

        //number of bits we use to encode pointers to the children
        size_t pt_bits = sym_width(node_n_bits/8);
        assert(pt_bits<(1<<phi_rep.int_pt_width));

        //amount of information (in bits) for the header
        size_t header_bits = 1+off_width+s_factor+phi_rep.int_pt_width+(pt_bits*n_children);
        header_bits = INT_CEIL(header_bits, 8)*8;//byte-aligned

        buffer.reserve_in_bits(header_bits+node_n_bits);

        //THE ENCODING STARTS HERE
        //==1 bit (false) to indicate this is an internal node
        size_t bit_pos = 0;
        buffer.write(bit_pos, bit_pos, 0);
        bit_pos++;
        //==

        //==offset information
        //rank_bits=r_width*node_sigma bits store the rank information
        buffer.write(bit_pos, bit_pos+off_width-1, node_offset);
        bit_pos+=off_width;
        //==

        //==
        //s_factor bits to indicate which children were collapsed
        //these bits give us the real number of children
        for(size_t i=0;i<s_factor;i++){
            buffer.write(bit_pos, bit_pos, child_marks[i]);
            bit_pos++;
        }
        assert(bit_pos==(1+off_width+s_factor));
        //==

        //==pt_width*n_children bits store pointers (byte offsets) to the children of this node
        buffer.write(bit_pos, bit_pos+phi_rep.int_pt_width-1, pt_bits);
        bit_pos+=phi_rep.int_pt_width;
        for(size_t i=0;i<n_children;i++){
            //NOTE consider the number of bytes in the header when computing the byte position of a child.
            //i.e., child_byte_pos = h_bits/8 + block_ptr[c], where c is the child we need to find
            buffer.write(bit_pos, bit_pos+pt_bits-1, block_ptr[i]);
            bit_pos+=pt_bits;
        }
        assert(bit_pos==(1+off_width+s_factor+phi_rep.int_pt_width+(pt_bits*n_children)));
        //==

        //move to the next byte-aligned position
        size_t header_bytes = INT_CEIL(bit_pos, 8);
        assert((header_bytes*8)==header_bits);

        //==add the stream of the children
        buffer.concatenate(header_bytes, children_buffer, node_n_bits/8);
        //==

        node_n_bits+=header_bits;

        //gather statistics
        stats.header_overhead+=header_bits;
    }

    inline void finish_tree() {

        assert(lvl==0);
        assert(aligned<8>(node_n_bits));

        tree_offset[n_children]= phi_rep.tot_syms;
        size_t n_blocks = INT_CEIL(phi_rep.tot_syms, b_size);//original number of blocks in the first level of the tree

        //pt_bits indicates how many bits we use to encode pointers to the trees:
        //node_n_bits/8 is the pointer and b_runs indicate collapsed blocks
        //so far, node_n_bits considers:
        // * the sum of the tree sizes in bits (excluding the external su/pred information)
        // * the sum of the ext. succ/pred information for the trees
        phi_rep.ext_pt_width = std::max<uint16_t>(sym_width(node_n_bits/8), 2*run_width)+1;

        //(pt_bits*n_blocks) for the pointers to the trees
        //+1 because we add a dummy tree at end for consistency
        size_t tree_ptr_bits = (phi_rep.ext_pt_width*n_blocks);
        size_t header_bits = tree_ptr_bits;
        phi_rep.header_bytes = INT_CEIL(header_bits, 8);
        header_bits = phi_rep.header_bytes*8;
        size_t bit_pos=header_bits;
        buffer.reserve_in_bits(header_bits+node_n_bits);

        auto *stream = (char *)buffer.stream;
        ifs->read(&stream[bit_pos/8], node_n_bits/8);

        //write the pointers to the trees
        bit_pos=0;
        size_t c=0, l=0, n_syms, r, offsets;
        for(size_t b=0;b<n_children;b++){

            buffer.write(bit_pos, bit_pos+phi_rep.ext_pt_width-1, (block_ptr[b]<<1));
            //std::cout<<"block:"<<b<<" real_block:"<<c<<" b_pos:"<<bit_pos<<" ptr:"<<block_ptr[b]<<" tree_offset:"<<tree_offset[b]<<" "<<tree_offset[b+1]<<std::endl;
            bit_pos+=phi_rep.ext_pt_width;
            r = (tree_offset[b+1]-tree_offset[b])/b_size;
            c++;
            l++;
            r--;
            n_syms =tree_offset[b];
            while((n_syms+b_size)<tree_offset[b+1]){
                offsets = (l<<run_width) | r;
                offsets = (offsets<<1) | 1;
                //std::cout<<"block:"<<b<<" real_block:"<<c<<" b_pos:"<<bit_pos<<" offsets:"<<l<<" "<<r<<std::endl;
                buffer.write(bit_pos, bit_pos+phi_rep.ext_pt_width-1, offsets);
                bit_pos+=phi_rep.ext_pt_width;
                n_syms +=b_size;
                c++;
                l++;
                r--;
            }
            l=0;
        }
        assert(bit_pos==tree_ptr_bits);
        assert(c==n_blocks);

        node_n_bits+=header_bits;

        //move the information to the bwt
        phi_rep.stream.swap(buffer);

        //gather statistics
        stats.header_overhead+=header_bits;
        stats.tree_pointers_overhead+=tree_ptr_bits;
    }

    void inline get_max_psum(uint64_t* tmp_psum, uint64_t* max_psum) const {
        //max prefix sum within 4 blocks with 8 runs each
        uint64_t tmp_max = std::max(std::max(tmp_psum[0], tmp_psum[1]),
                                    std::max(tmp_psum[2], tmp_psum[3]));
        if(tmp_max>max_psum[0]){
            max_psum[0] = tmp_max;
        }

        //max prefix sum within 2 blocks with 16 runs each
        tmp_psum[0]+=tmp_psum[1];
        tmp_psum[2]+=tmp_psum[3];
        tmp_max = std::max(tmp_psum[0], tmp_psum[2]);
        if(tmp_max>max_psum[1]){
            max_psum[1] = tmp_max;
        }

        //max prefix sum within one blocks withn 32 runs
        tmp_psum[0]+=tmp_psum[2];
        if(tmp_max>max_psum[2]){
            max_psum[2] = tmp_psum[0];
        }
    }

    inline uint8_t compute_leaf_enc_code(uint8_t max_bytes, bool vbyte_enc, const uint64_t *max_psum) const {

        uint8_t code = 0;

        switch (max_bytes) {
            case 1://1-byte encoding (we can read 16 (or 32) runs at the time)
                code +=max_psum[1]>0xFF;//overflow with 16 runs?
                code +=max_psum[2]>0XFF;//overflow with 32 runs?
                return code;
            case 2://2-byte encoding (we can read 8 (or 16) runs at the time)
                code =3;
                code +=max_psum[0]>0xFFFF;//overflow with 8 runs?
                code +=max_psum[1]>0xFFFF;//overflow with 16 runs?
                code <<=vbyte_enc;
                return code;
            case 3://3-byte encoding (we can read 4 (or 8) runs at the time. We assume no overflow)
                code =9;
                code +=vbyte_enc;//using vbytes?
                return code;
            case 4://four-byte encoding (we can read 4 (or 8) runs at the time. We assume no overflow)
                code =11;
                code +=vbyte_enc;//using vbytes?
                return code;
            case 5://five-byte encoding (we can read 2 (or 4) runs at the time. We assume no overflow
                return 13;
            case 6:
                return 14;
            case 7:
                return 15;
            default:
                std::cout<<"Unknown code for leaf encoding"<<std::endl;
                exit(1);
        }
    }

    inline void create_leaf(std::vector<block_type>& blocks, size_t n_blocks){

        //for 1 byte: check that the sum of 16 (or 32 for AVX) consecutive run lens is <=256
        //for 2 bytes: check that the sum of 8 (or 16 for AVX) consecutive run lens is <=2^16-1
        //for 3-4 bytes: check that the sum of 4 (or 8 for AVX) consecutive run lens is <=2^32-1
        //for 4-8 bytes: check that the sum of 2 (or 4 for AVX) consecutive run lens is <=2^64-1
        //we need to do this check to work with SIMD instructions. If it does not fit, use the next encoding that can fits them

        assert(node_n_bits==0);
        assert(lvl>0);

        size_t sym, len, n_runs=0;
        n_runs+=blocks[0].size()-1;
        //combine the runs we previously broke due to boundary constraints
        for(size_t i=1;i<n_blocks;i++){
            //read the rightmost run of the previous block
            sym = blocks[i-1].back().first;
            len = blocks[i-1].back().second;

            //collapse the run with the first run of the current block if they have the same symbol
            if(sym==blocks[i][0].first){
                //collapse the runs as they are the same
                blocks[i][0].second +=len;
                blocks[i-1].pop_back();
            } else {
                //otherwise process the last run of the previous block as an independent run
                n_runs++;
            }
            n_runs+=blocks[i].size()-1;
        }
        //the last run
        n_runs++;

        //compute the total number of vbytes the run lengths use and
        // if there are prefix sums causing overflow
        assert(n_runs<=phi_dt_type::max_block_runs);
        size_t bfr_dist[9]={0}, bytes;
        uint64_t tmp_psum[4]={0};//8,16,24,32
        uint64_t max_psum[3]={0};//8,16,32
        uint64_t bk = 0;
        uint64_t max_sym=0;
        for(size_t i=0;i<n_blocks;i++){
            for(auto & run : blocks[i]){
                // get the frequency of bytes that the run lengths use
                // (they should be relatively small numbers. Most of them should fit 1-2 bytes)
                bytes = INT_CEIL(sym_width(run.second), 8);
                bfr_dist[bytes]++;
                tmp_psum[bk>>3] += run.second;
                bk++;

                if(bk==32){
                    get_max_psum(tmp_psum, max_psum);
                    memset(tmp_psum, 0, 32);
                    bk=0;
                }
                if(run.first>max_sym) max_sym = run.first;
            }
        }
        get_max_psum(tmp_psum, max_psum);
        assert(bfr_dist[0]==0);

        size_t total_vbytes=0, max_bytes=0;
        for(size_t b=1;b<9;b++){
            total_vbytes +=bfr_dist[b]*b;
            if(bfr_dist[b]>0) max_bytes = b;
        }
        assert(max_bytes>0 && max_bytes<6);
        //

        if(max_bytes>1){
            //number of control masks of 1 byte for fast vbyte decoding;
            total_vbytes += INT_CEIL(n_runs, (8/sym_width(max_bytes-1)));
        }

        //alternative encoding using a fixed number of bytes per run
        size_t total_fbytes = max_bytes*n_runs;

        //choose the byte encoding for the run lengths of this leaf
        size_t len_bits;
        bool fix_len_enc=false;
        if(total_vbytes<total_fbytes){
            assert(max_bytes>1);
            len_bits = total_vbytes*8;
        } else {
            len_bits = total_fbytes*8;
            fix_len_enc = true;
        }
        size_t sym_bits = sym_width(max_sym)*n_runs;
        sym_bits = INT_CEIL(sym_bits, 8)*8;

        //
        len_enc = compute_leaf_enc_code(max_bytes, !fix_len_enc, max_psum);

        //==HEADER DESCRIPTION:
        //1 bit to indicate it is a leaf node
        //phi_rep.len_enc_width bits to store the encoding of the run lengths
        //32 bits to store the value of (run_bits/8)
        //phi_rep.int_pt_width bits to store the value of sym_width(max_sym)
        //===

        size_t header_bits = 1+phi_rep.len_enc_width+32+phi_rep.int_pt_width;
        header_bits = INT_CEIL(header_bits, 8)*8;//byte-aligned

        //allocate bytes for the information of this leaf
        buffer.reserve_in_bits(header_bits + len_bits + sym_bits);

        //start writing the in the buffer
        size_t bit_pos = 0;
        //1 bit (true) to indicate this node is a leaf
        buffer.write(bit_pos, bit_pos, 1);
        bit_pos++;

        //store the encoding of the run lengths
        buffer.write(bit_pos, bit_pos+phi_rep.len_enc_width-1, len_enc);
        bit_pos+=phi_rep.len_enc_width;
        //

        //store the number of bytes used by the run lengths
        buffer.write(bit_pos, bit_pos+31, len_bits/8);
        bit_pos+=32;
        //

        //store the number of bits we use to encode the symbols
        buffer.write(bit_pos, bit_pos+phi_rep.int_pt_width-1, sym_width(max_sym));
        bit_pos+=phi_rep.int_pt_width;

        //the run lengths are byte-aligned
        size_t byte_pos = INT_CEIL(bit_pos, 8);
        assert((byte_pos*8)==header_bits);

        auto *byte_stream = (uint8_t *) buffer.stream;
        size_t written_bytes = insert_run_lens(blocks, n_blocks, &byte_stream[byte_pos], max_bytes , fix_len_enc);
        assert((written_bytes*8)==len_bits);

        bit_pos = (byte_pos*8)+len_bits;
        uint8_t w = sym_width(max_sym);
        for(size_t i=0;i<n_blocks;i++) {
            for (auto &run: blocks[i]) {
                buffer.write(bit_pos, bit_pos+w-1, run.first);
                bit_pos+=w;
            }
        }
        bit_pos = INT_CEIL(bit_pos, 8)*8;
        assert(bit_pos==(header_bits+len_bits+sym_bits));

        node_n_bits = header_bits+len_bits+sym_bits;
        phi_rep.eff_runs += n_runs;

        //gather statistics
        if(n_blocks>stats.max_n_blocks) stats.max_n_blocks = n_blocks;

        stats.sym_overhead+=sym_bits;
        stats.len_overhead+=len_bits;
        stats.runs_overhead+=len_bits+sym_bits;
        stats.header_overhead+=header_bits;
        stats.rpl_freq[n_runs]++;
        stats.leaf_depth_freq[lvl-1]++;//lvl=0 is the tree, so it doesn't count. lvl=1 is a root of a block
        stats.leaf_enc_freq[len_enc]++;//the encoding type for a leaf
    }

    size_t insert_run_lens(std::vector<block_type>& blocks, size_t n_blocks, uint8_t *stream, size_t max_bytes, bool fix_len_enc){

        size_t written_bytes=0;
        if(fix_len_enc){
            size_t enc_run;
            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){
                    enc_run =  run.second;
                    memcpy(stream, &enc_run, max_bytes);
                    stream+=max_bytes;
                    written_bytes+=max_bytes;
                }
            }
        }else{

            uint8_t n_ctrl_bits = sym_width(max_bytes-1);
            size_t vb_lens[8]={0};

            //we pack the stream in groups of (at most) 8 elements,
            // Each code uses (at most) 8 bytes, thus we need 8*8=64 tmp_bytes
            uint64_t code;
            uint8_t tmp_stream[64];
            uint8_t ctrl_bits = 0, acc_width=0, p=0, byte_pos=0;

            for(size_t i=0;i<n_blocks;i++){
                for(auto & run : blocks[i]){

                    vb_lens[p] = INT_CEIL(sym_width(run.second), 8);
                    code = run.second;
                    memcpy(&tmp_stream[byte_pos], &code, vb_lens[p]);
                    ctrl_bits |= ((vb_lens[p]-1U) << acc_width);
                    byte_pos+=vb_lens[p];
                    p++;
                    acc_width+=n_ctrl_bits;

                    if(acc_width+n_ctrl_bits>8){
                        *stream=ctrl_bits;
                        stream++;

                        memcpy(stream, &tmp_stream[0], byte_pos);
                        stream+=byte_pos;
                        written_bytes+=byte_pos+1;
                        p=0;
                        ctrl_bits = 0;
                        acc_width = 0;
                        byte_pos = 0;
                    }
                }
            }

            if(p!=0){
                *stream=ctrl_bits;
                stream++;
                memcpy(stream, &tmp_stream[0], byte_pos);
                written_bytes+=byte_pos+1;
            }
        }
        return written_bytes;
    }

    void print_node_info(std::vector<block_type> bkl, size_t n_blocks){

        std::string pad = std::string(lvl+1, '\t');

        std::cout<<pad<<(leaf? "leaf" :"internal_node")<<std::endl;
        std::cout<<pad<<"node_offset:"<<node_offset<<std::endl;
        std::cout<<pad<<"lm_child:"<<lm_child<<std::endl;
        std::cout<<pad<<"rm_child:"<<rm_child<<std::endl;
        std::cout<<pad<<"sym_before:"<<syms_before<<std::endl;
        std::cout<<pad<<"covered_symbols:"<<cov_symbols<<std::endl;
        std::cout<<pad<<"child_rank:"<<child_rank<<std::endl;
        std::cout<<pad<<"level:"<<int(lvl)<<std::endl;
        std::cout<<pad<<"size in bytes:"<<INT_CEIL(node_n_bits, 8)<<std::endl;
        std::cout<<pad<<"lm_branch:"<<lm_tree_branch<<std::endl;
        std::cout<<pad<<"rm_branch:"<<rm_tree_branch<<std::endl;

        if(!leaf){
            std::cout<<pad<<"n_children:"<<n_children<<std::endl;
            std::cout<<pad<<"child info:";
            for(size_t k=0;k<s_factor;k++){
                std::cout<<(k>0 ? ", ":"")<<int(child_marks[k]);
            }
            std::cout<<""<<std::endl;
            std::cout<<pad<<"pointer to children:";
            for(size_t k=0;k<n_children;k++){
                std::cout<<(k>0 ? ", ":"")<<block_ptr[k];
            }
        }else{
            std::cout<<pad<<"leaf encoding:"<<int(len_enc)<<std::endl;
            std::cout<<pad<<"runs: ";
            for(size_t k=0;k<n_blocks;k++){
                for(auto & r : bkl[k]){
                    std::cout<<"(sym:"<<r.first<<",len:"<<r.second<<") ";
                }
            }
        }
        std::cout<<""<<std::endl;
    }

    template<node_type type>//internal or leaf
    inline void create_node(size_t n_blocks) {

        assert(aligned<8>(node_n_bits));//check it is byte-aligned
        //lvl=0 means the tree root, and tmp_node is then the root v of a block in the tree.
        // Therefore, v is the leftmost and rightmost branches of the tree
        tmp_node->cov_symbols = b_size*n_blocks;
        tmp_node->lm_child = lvl==0 || consumed_syms == 0;
        tmp_node->rm_child = lvl==0 || (consumed_syms+ tmp_node->cov_symbols)==(b_size*s_factor);
        tmp_node->lm_tree_branch = lm_tree_branch && tmp_node->lm_child;
        tmp_node->rm_tree_branch = rm_tree_branch && tmp_node->rm_child;
        tmp_node->child_rank = n_children;
        tmp_node->leaf = type==LEAF;
        tmp_node->syms_before=syms_before+consumed_syms;
        tmp_node->node_offset = consumed_syms;

        if constexpr (type==INTERNAL){
            assert(n_blocks==1);
            for(auto & run : active_blocks[0]){
                tmp_node->process_run(run.first, run.second);
            }
            tmp_node->finish_run_scan();
            tmp_node->finish_int_node();
            stats.children_freq[tmp_node->n_children]++;
        } else {
            assert(n_blocks>=1);
            tmp_node->create_leaf(active_blocks, n_blocks);
        }

        //print the node information for debugging purposes
        //tmp_node->print_node_info(active_blocks, n_blocks);
        //
        assert(aligned<8>(node_n_bits+tmp_node->node_n_bits));

        if(lvl>0){
            child_marks[child_mark_acc]=true;
            child_mark_acc+=n_blocks;
            assert(child_mark_acc<=s_factor);
            //append the stream of tmp_node to the stream of children for this node
            children_buffer.reserve_in_bits(node_n_bits+tmp_node->node_n_bits);
            children_buffer.concatenate(node_n_bits/8, tmp_node->buffer, tmp_node->node_n_bits/8);
        } else {
            tree_offset[n_children] = consumed_syms;
            //store to disk
            assert(ofs->tellp()==block_ptr[n_children]);
            assert(aligned<8>(tmp_node->node_n_bits));
            ofs->write((char *)tmp_node->buffer.stream, tmp_node->node_n_bits/8);
            stats.trees_overhead+=tmp_node->node_n_bits;
        }

        //the pointer to the active child node (next_node) should be aligned
        node_n_bits+=tmp_node->node_n_bits;
        block_ptr[++n_children]=(node_n_bits/8);
        consumed_syms+=tmp_node->cov_symbols;
        tmp_node->reset();
    }

    inline void reset(){
        std::fill(child_marks.begin(), child_marks.end(), false);
        n_children = 0;
        node_n_bits = 0;
        consumed_syms = 0;
        lm_tree_branch = lvl==0;//lvl=1 means the root of the tree
        lm_tree_branch = lvl==0;
        block_ptr[0] = 0;
        child_mark_acc=0;
        len_enc=0;
        syms_before=0;
    }
};

template<class phi_dt_type, class size_type>
struct phi_tree{

    typedef phi_node<phi_dt_type, size_type> node_type;
    typedef typename node_type::run_type run_type;

    phi_stat_collector<phi_dt_type> stats;
    tmp_workspace twd;
    phi_dt_type& phi_rep;
    node_type *root= nullptr;

    explicit phi_tree(const std::string& tmp_dir, phi_dt_type& _phi_rep): twd(tmp_dir),
                                                                          phi_rep(_phi_rep){}

    void build(std::string& rsa_samp_file) {
        //compute basic statistics
        size_t f_size = std::filesystem::file_size(rsa_samp_file);
        size_t n_runs = f_size / sizeof(run_type);
        size_t buff_size = 1024*1024;
        size_t n_blocks = n_runs/buff_size;
        size_t rem_samples = n_runs;
        std::vector<run_type> rsa_samples_buff(buff_size);
        std::ifstream rsa_samp_ifs(rsa_samp_file, std::ios::binary);
        for(size_t i=0;i<n_blocks;i++){
            rsa_samp_ifs.read((char *)rsa_samples_buff.data(), buff_size*sizeof(run_type));
            for(size_t j=0;j<buff_size;j++){
                phi_rep.tot_syms+=rsa_samples_buff[j].second;
            }
            rem_samples=-buff_size;
        }
        rsa_samp_ifs.read((char *)rsa_samples_buff.data(), rem_samples*sizeof(run_type));
        for(size_t j=0;j<rem_samples;j++){
            phi_rep.tot_syms+=rsa_samples_buff[j].second;
        }
        phi_rep.orig_runs = n_runs;
        //

        root = new node_type(0, phi_dt_type::block_size, phi_rep, stats);
        std::ofstream ofs(twd.get_file("trees"), std::ios::binary);
        root->ofs = &ofs;

        node_type *current = root;
        size_t b_size= phi_dt_type::block_size/phi_dt_type::scale_factor;
        for(size_t i=1;i<=phi_rep.levels;i++){
            current->tmp_node = new node_type(i, b_size, phi_rep, stats);
            current = current->tmp_node;
            b_size/=phi_dt_type::scale_factor;
        }

        //compute the tree
        for(size_t i=0;i<n_blocks;i++){
            rsa_samp_ifs.read((char *)rsa_samples_buff.data(), buff_size*sizeof(run_type));
            for(size_t j=0;j<buff_size;j++){
                root->process_run(rsa_samples_buff[j].first, rsa_samples_buff[j].second);
            }
            rem_samples=-buff_size;
        }
        rsa_samp_ifs.read((char *)rsa_samples_buff.data(), rem_samples*sizeof(run_type));
        for(size_t j=0;j<rem_samples;j++){
            root->process_run(rsa_samples_buff[j].first, rsa_samples_buff[j].second);
        }
        rsa_samp_ifs.close();
        destroy_vector(rsa_samples_buff);

        root->finish_run_scan();
        ofs.close();

        std::ifstream ifs(twd.get_file("trees"), std::ios::binary);
        root->ifs = &ifs;
        root->finish_tree();
        ifs.close();
        //
    }

    void build_in_memory(std::vector<run_type>& block){

        //basic information about the block (the area it covers (number of symbols) and the number of runs)
        phi_rep.orig_runs = block.size();
        phi_rep.tot_syms=0;
        for(size_t j=0;j<block.size();j++){
            phi_rep.tot_syms += block[j].second;
        }
        //

        root = new node_type(0, phi_dt_type::block_size, phi_rep, stats);
        std::ofstream ofs(twd.get_file("trees"), std::ios::binary);
        root->ofs = &ofs;

        node_type *current = root;
        size_t b_size= phi_dt_type::block_size/phi_dt_type::scale_factor;
        for(size_t i=1;i<=phi_rep.levels;i++){
            current->tmp_node = new node_type(i, b_size, phi_rep, stats);
            current = current->tmp_node;
            b_size/=phi_dt_type::scale_factor;
        }

        //compute the tree
        for(size_t j=0;j<block.size();j++){
            root->process_run(block[j].first, block[j].second);
        }
        root->finish_run_scan();
        ofs.close();

        std::ifstream ifs(twd.get_file("trees"), std::ios::binary);
        root->ifs = &ifs;
        root->finish_tree();
        ifs.close();
        //
    }

    void report_stats(){

        std::cout<<"Number_of_runs_in_a_leaf dist:"<<std::endl;
        assert(stats.rpl_freq[0]==0);
        for(size_t r=1;r<=phi_dt_type::max_block_runs;r++){
            std::cout<<"\t"<<r<<" : "<<stats.rpl_freq[r]<<std::endl;
        }
        std::cout<<"Total number of runs versus original number of runs: "<<phi_rep.eff_runs<<" / "<<phi_rep.orig_runs<<std::endl;
        std::cout<<"Increase in the number of runs "<<(double(phi_rep.eff_runs)/double(phi_rep.orig_runs)-1)*100<<"%"<<std::endl;

        std::cout<<"Leaf_depth dist:"<<std::endl;
        size_t tot_leaves=0;
        for(size_t i=0;i<20;i++){
            tot_leaves+=stats.leaf_depth_freq[i];
        }
        for(size_t i=0;i<20;i++){
            if(stats.leaf_depth_freq[i]>0){
                std::cout<<"\t"<<i<<": "<<double(stats.leaf_depth_freq[i])/double(tot_leaves)<<std::endl;
            }
        }

        std::cout<<"Leaf_encoding dist:"<<std::endl;
        for(size_t i=0;i<16;i++){
            switch(i) {
                case 0:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 1:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte, overflow 16: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 2:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t1 byte, overflow 32: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 3:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 4:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 8: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 5:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 16: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 6:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 7:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 8, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 8:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t2 bytes, overflow 16, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 9:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t3 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 10:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t3 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 11:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t4 bytes: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 12:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t4 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 13:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t5 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                case 14:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t6 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
                    break;
                default:
                    if(stats.leaf_enc_freq[i]) std::cout<<"\t7 bytes, vbyte_comp: "<<double(stats.leaf_enc_freq[i])/double(tot_leaves)<<std::endl;
            }
        }

        std::cout<<"Number_of_children dist:"<<std::endl;
        size_t del_nodes=0, tot_nodes=0;
        for(size_t i=0;i<20;i++){
            if(stats.children_freq[i]!=0){
                std::cout<<"\t"<<i<<": "<<stats.children_freq[i]<<std::endl;
                del_nodes+=(phi_rep.scale_factor-i)*stats.children_freq[i];
                tot_nodes+=stats.children_freq[i];
            }
        }

        size_t n_blocks = INT_CEIL(phi_rep.tot_syms, phi_rep.block_size);//original number of blocks in the first level of the tree
        std::cout<<"Effective number of trees versus full number of trees (n/b): "<<root->n_children<<" / "<<n_blocks<<std::endl;
        std::cout<<"Percentage of removed trees: "<<(1-double(root->n_children)/double(n_blocks))*100<<"% "<<std::endl;

        tot_nodes*=phi_rep.scale_factor;
        tot_nodes+=n_blocks;
        del_nodes=n_blocks-root->n_children;

        std::cout<<"Percentage of removed nodes: "<<(double(del_nodes)/double(tot_nodes))*100<<"% "<<std::endl;
        std::cout<<"Written bytes in the data structure: "<<INT_CEIL(root->node_n_bits, 8)<<std::endl;
        std::cout<<"Space breakdown"<<std::endl;
        std::cout<<"\tRuns: "<<INT_CEIL(stats.runs_overhead, 8)<<" bytes ("<<(double(stats.runs_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tSymbols overhead: "<<INT_CEIL(stats.sym_overhead, 8)<<" bytes ("<<(double(stats.sym_overhead)/double(stats.runs_overhead))*100<<"%)"<<std::endl;
        std::cout<<"\t\tLengths overhead: "<<INT_CEIL(stats.len_overhead, 8)<<" bytes ("<<(double(stats.len_overhead)/double(stats.runs_overhead))*100<<"%)"<<std::endl;
        std::cout<<"\tHeaders: "<<INT_CEIL(stats.header_overhead, 8)<<" bytes ("<<(double(stats.header_overhead)/double(root->node_n_bits))*100<<"%)"<<std::endl;
        std::cout<<"\t\tTree pointers: "<<INT_CEIL(stats.tree_pointers_overhead, 8)<<" bytes ("<<(double(stats.tree_pointers_overhead)/double(stats.header_overhead))*100<<"%)"<<std::endl;
        size_t ptr_bv_ov = stats.header_overhead - stats.tree_pointers_overhead;
        std::cout<<"\t\tInt. Pointers, bitvectors, and extras: "<<INT_CEIL(ptr_bv_ov, 8)<<" bytes ("<<(double(ptr_bv_ov)/double(stats.header_overhead))*100<<"%)"<<std::endl;
        assert((stats.header_overhead+stats.runs_overhead)==root->node_n_bits);
        //std::cout<<"\t\tTrees pointers: "<<INT_CEIL(stats.trees_overhead, 8)<<std::endl;
        std::cout<<"space_usage:"<<float(root->node_n_bits)/float(phi_rep.tot_syms)<<" bps"<<std::endl;
    }
};

template<class phi_dt_type, class size_type, bool vbyte=false>
void build_vlbt_phi(phi_dt_type& phi_rep, std::string& rsa_samp_file, std::string tmp_dir="./"){
    phi_tree<phi_dt_type, size_type> tree(tmp_dir, phi_rep);
    tree.build(rsa_samp_file);
    tree.report_stats();
}

template<class phi_dt_type, class size_type, bool vbyte=false>
void build_vlbt_phi_in_memory(phi_dt_type& phi_rep,
                              std::vector<std::pair<size_type, size_type>> block,
                              std::string tmp_dir="./"){
    phi_tree<phi_dt_type, size_type> tree(tmp_dir, phi_rep);
    tree.build_in_memory(block);
    tree.report_stats();
}

struct rsa_type {
    uint64_t tail_val;
    uint64_t next_head_val;
    uint64_t run;
};

uint64_t get_diff(uint64_t first, uint64_t second) {
    uint64_t abs_diff = (first > second) ? (first - second): (second - first);
    assert(abs_diff<=INT64_MAX);
    return abs_diff;
}

template<class size_type>
void preprocess_rsa(std::string& rsa_file, std::string& rsa_per_str_file,
                    size_t ssamp_val, std::string& ssamp_phi_file, std::string& ssamp_th_file){

    size_t n_elements = std::filesystem::file_size(rsa_file)/sizeof(size_type);
    std::cout<<"There are "<<n_elements<<" SA samples"<<std::endl;
    std::vector<rsa_type> samples(n_elements/2);
    size_t s_pos=0;

    std::ifstream ifs(rsa_file, std::ios::binary);
    size_t buffer_size = 4096;
    std::vector<size_type> buffer(buffer_size+1, 0);
    size_t n_blocks = n_elements/buffer_size;
    size_t rem = n_elements;

    //read samples from disk and reorganize them
    ifs.read((char *)buffer.data(), off_t(sizeof(size_type)*buffer_size));
    for(size_t i=0;i<(n_blocks-1);i++){
        for(size_t j=1;j<buffer_size;j+=2){
            samples[s_pos].tail_val = buffer[j];
            samples[s_pos].next_head_val = buffer[j+1];
            samples[s_pos].run = s_pos++;
        }
        assert(samples[s_pos-1].next_head_val==0);
        rem -=buffer_size;
        ifs.read((char *)buffer.data(), off_t(sizeof(size_type)*buffer_size));
        samples[s_pos-1].next_head_val = buffer[0];
    }

    rem -=buffer_size;
    for(size_t j=1;j<buffer_size;j+=2){
        samples[s_pos].tail_val = buffer[j];
        samples[s_pos].next_head_val = buffer[j+1];
        samples[s_pos].run = s_pos++;
    }
    assert(samples[s_pos-1].next_head_val==0);

    if(rem>0){
        ifs.read((char *)buffer.data(), off_t(sizeof(size_type)*rem));
        for(size_t j=1;j<rem;j+=2){
            samples[s_pos].tail_val = buffer[j];
            samples[s_pos].next_head_val = buffer[j+1];
            samples[s_pos].run = s_pos++;
        }
    }
    ifs.close();
    assert(s_pos==samples.size());
    assert(s_pos==(n_elements/2));
    //

    //sort the samples by text position
    uint64_t last_run = s_pos-1;
    std::sort(samples.begin(), samples.end(), [](auto const& a, auto const&b){
        return a.tail_val<b.tail_val;
    });
    //

    //compute and store the subsamples for phi
    size_t f_size = std::filesystem::file_size(rsa_per_str_file);
    n_elements = f_size/sizeof(size_type);
    std::vector<size_type> str_ranges(n_elements, 0);
    std::ifstream ifs2(rsa_per_str_file, std::ios::binary);
    ifs2.read((char *)str_ranges.data(), off_t(f_size));
    uint64_t discard_mark = std::numeric_limits<uint64_t>::max();
    std::ofstream ifs_phi(ssamp_phi_file, std::ios::binary);

    size_t best_comp=0;

    size_t n_strings = n_elements-1, last_sampled, n_samp=0, len, acc_len=0, buff_pos=0;
    size_type str_boundary;
    s_pos=0;
    for(size_t str=0;str<n_strings;str++){
        str_boundary = str_ranges[str+1]-1;
        last_sampled = s_pos;
        n_samp++;
        s_pos++;
        while(samples[s_pos].tail_val<str_boundary){
            if((samples[s_pos+1].tail_val-samples[last_sampled].tail_val>ssamp_val) ||
                samples[s_pos].run==last_run){

                len = samples[s_pos].tail_val-samples[last_sampled].tail_val;
                buffer[buff_pos++] = samples[last_sampled].next_head_val;
                buffer[buff_pos++] = len;
                best_comp+=sym_width(get_diff(samples[last_sampled].next_head_val, samples[last_sampled].tail_val));

                if(buff_pos==buffer_size){
                    ifs_phi.write((char *)buffer.data(), sizeof(size_type)*buffer_size);
                    buff_pos=0;
                }
                acc_len+=len;
                //std::cout<<"tail_pos:"<<samples[last_sampled].tail_val<<", next_head_val:"<<samples[last_sampled].next_head_val<<", run:"<<samples[last_sampled].run<<std::endl;
                //std::cout<<"run:("<<samples[last_sampled].next_head_val<<","<<len<<")"<<std::endl;
                last_sampled = s_pos;
                n_samp++;
            } else {
                samples[s_pos].next_head_val = discard_mark;
            }
            s_pos++;
        }
        len = (str_boundary+1)-samples[last_sampled].tail_val;
        buffer[buff_pos++] = samples[last_sampled].next_head_val;
        buffer[buff_pos++] = len;
        best_comp+=sym_width(get_diff(samples[last_sampled].next_head_val, samples[last_sampled].tail_val));
        if(buff_pos==buffer_size){
            ifs_phi.write((char *)buffer.data(), sizeof(size_type)*buffer_size);
            buff_pos=0;
        }
        acc_len += len;
        //std::cout<<"tail_pos:"<<samples[last_sampled].tail_val<<", next_head_val:"<<samples[last_sampled].next_head_val<<", run:"<<samples[last_sampled].run<<std::endl;
        //std::cout<<"run:("<<samples[last_sampled].next_head_val<<","<<len<<")"<<std::endl;
    }
    //
    std::cout<<"The best compression we can achieve is "<<INT_CEIL(best_comp, 8)<<" bytes "<<std::endl;
    assert(acc_len==str_ranges.back());
    if(buff_pos>0){
        ifs_phi.write((char *)buffer.data(), sizeof(size_type)*buff_pos);
    }
    ifs_phi.close();
    //store them as runs
    /*n_samp=0; s_pos=0;
    for(size_t str=0;str<n_strings;str++){
        str_boundary = str_ranges[str+1]-1;
        last_sampled = s_pos;
        n_samp++;
        s_pos++;
        while(samples[s_pos].tail_val<str_boundary){
            if((samples[s_pos+1].tail_val-samples[last_sampled].tail_val>sub_samp_val) ||
               samples[s_pos].run==last_run){
                last_sampled = s_pos;
                n_samp++;
            }else{
                samples[s_pos].next_head_val = discard_mark;
            }
            s_pos++;
        }
    }*/
    //
    /*std::cout<<" ==== sampled ==== "<<std::endl;
    for(size_t i=0;i<20;i++){
        if(samples[i].next_head_val!=discard_mark){
            std::cout<<samples[i].tail_val<<" "<<samples[i].next_head_val<<std::endl;
        }
    }*/
    std::cout<<"We subsampled "<<n_samp<<" elements out of "<<samples.size()<<" ("<<double(n_samp)/double(samples.size())*100<<"%)"<<std::endl;
    //
}

#endif//VLBT_CONSTRUCT_VLBT_PHI