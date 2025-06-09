
// Extracts random patterns from a file

/*#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

static int Seed;
#define ACMa 16807
#define ACMm 2147483647
#define ACMq 127773
#define ACMr 2836
#define hi (Seed / ACMq)
#define lo (Seed % ACMq)

static int fst = 1;

/ *
* returns a random integer in 0..top-1
* /
int aleat(int top) {
	long test;
	struct timeval t;
	if (fst) {
		gettimeofday (&t, NULL);
		Seed = t.tv_sec * t.tv_usec;
		fst = 0;
	}
	{
		Seed = ((test = ACMa * lo - ACMr * hi) > 0) ? test : test + ACMm;
		return ((double) Seed) * top / ACMm;
	}
}

void parse_forbid(unsigned char *forbid, unsigned char ** forbide) {

	int len, i, j;
	len = strlen(forbid);
	
	*forbide = (unsigned char *) malloc((len+1)*sizeof(unsigned char));
	if (*forbide == NULL)
	{
		fprintf (stderr, "Error: cannot allocate %i bytes\n", len+1);
		fprintf (stderr, "errno = %i\n", errno);
		exit (1);
	}

	for(i = 0, j = 0; i < len; i++) {
		if(forbid[i] != '\\') {
			if(forbid[i] != '\n')
				(*forbide)[j++] = forbid[i];
		} else { 
			i++;
			if(i == len) {
				forbid[i-1] = '\0';
				(*forbide)[j] = '\0';
				fprintf (stderr, "Not correct forbidden string: only one \\\n");
				return;
			}
			switch (forbid[i]) {
				case'n':  (*forbide)[j++] = '\n'; break;
				case'\\': (*forbide)[j++] = '\\'; break;
				case'b':  (*forbide)[j++] = '\b'; break;				
				case'e':  (*forbide)[j++] = '\e'; break;
				case'f':  (*forbide)[j++] = '\f'; break;
				case'r':  (*forbide)[j++] = '\r'; break;
				case't':  (*forbide)[j++] = '\t'; break;
				case'v':  (*forbide)[j++] = '\v'; break;
				case'a':  (*forbide)[j++] = '\a'; break;
				case'c':  
						if(i+3 >= len) {
							forbid[i-1] = '\0';
							(*forbide)[j] = '\0';
							fprintf (stderr, "Not correct forbidden string: 3 digits after \\c\n");
							return;
						}
						(*forbide)[j++] = (forbid[i+1]-48)*100 +
										  (forbid[i+2]-48)*10 + (forbid[i+3]-48); 
						i+=3;
						break;					
				default:
				fprintf (stdout, "Unknown escape sequence '\\%c'in forbidden string\n", forbid[i]);
				break;
			}
		}
	}
	(*forbide)[j] = '\0';
}

int main (int argc, char **argv) {
	int n, m, J, t;
	struct stat sdata;
	FILE *ifile, *ofile;
	unsigned char *buff;
	unsigned char *forbid, *forbide = NULL;

	if (argc < 5) {
		fprintf (stderr,
			 "Usage: genpatterns <file> <length> <number> <patterns file> <forbidden>\n"
			 "  randomly extracts <number> substrings of length <length> from <file>,\n"
			 "  avoiding substrings containing characters in <forbidden>.\n"
			 "  The output file, <patterns file> has a first line of the form:\n"
			 "    # number=<number> length=<length> file=<file> forbidden=<forbidden>\n"
			 "  and then the <number> patterns come successively without any separator.\n"
			 "  <forbidden> uses \\n, \\t, etc. for nonprintable chracters or \\cC\n"
			 "  where C is the ASCII code of the character written using 3 digits.\n\n");
		exit (1);
	}

	if (stat (argv[1], &sdata) != 0) {
		fprintf (stderr, "Error: cannot stat file %s\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}
	n = sdata.st_size;

	m = atoi (argv[2]);
	if ((m <= 0) || (m > n)) {
		fprintf (stderr,
			 "Error: length must be >= 1 and <= file length"
			 " (%i)\n", n);
		exit (1);
	}

	J = atoi (argv[3]);
	if (J < 1) {
		fprintf (stderr, "Error: number of patterns must be >= 1\n");
		exit (1);
	}

	if (argc > 5) {
		forbid = argv[5];
		parse_forbid(forbid, &forbide);
	} else
		forbid = NULL;

	ifile = fopen (argv[1], "r");
	if (ifile == NULL) {
		fprintf (stderr, "Error: cannot open file %s for reading\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	buff = (unsigned char *) malloc (n);
	if (buff == NULL) {
		fprintf (stderr, "Error: cannot allocate %i bytes\n", n);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	if (fread (buff, n, 1, ifile) != 1) {
		fprintf (stderr, "Error: cannot read file %s\n", argv[1]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}
	fclose (ifile);

	ofile = fopen (argv[4], "w");
	if (ofile == NULL) {
		fprintf (stderr, "Error: cannot open file %s for writing\n",
			 argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	if (fprintf (ofile, "# number=%i length=%i file=%s forbidden=%s\n",
		     J, m, argv[1],
		     forbid == NULL ? "" : (char *) forbid) <= 0) {
		fprintf (stderr, "Error: cannot write file %s\n", argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	for (t = 0; t < J; t++) {
		int j, l;
		if (!forbide)
			j = aleat (n - m + 1);
		else
		{
			do
			{
				j = aleat (n - m + 1);
				for (l = 0; l < m; l++)
					if (strchr (forbide, buff[j + l]))
						break;
			}
			while (l < m);
		}
		for (l = 0; l < m; l++)
			if (putc (buff[j + l], ofile) != buff[j + l])
			{
				fprintf (stderr,
					 "Error: cannot write file %s\n",
					 argv[4]);
				fprintf (stderr, " errno = %i\n", errno);
				exit (1);
			}
	}

	if (fclose (ofile) != 0) {
		fprintf (stderr, "Error: cannot write file %s\n", argv[4]);
		fprintf (stderr, " errno = %i\n", errno);
		exit (1);
	}

	fprintf (stderr, "File %s successfully generated\n", argv[4]);
	free(forbide);
    return 0;
}*/

#include <cassert>
#include <string>
#include <iostream>
#include <random>
#include <unordered_set>
#include <filesystem>
#include <fstream>

std::vector<uint64_t> sample_random_positions(std::unordered_set<uint64_t>& seen_positions, uint64_t n, uint64_t x) {

    if (x > n) throw std::invalid_argument("x cannot be larger than n");

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, n - 1);
    std::unordered_set<uint64_t> new_positions;

    while(new_positions.size() < x) {
        uint64_t val = dist(rng);
        if(seen_positions.find(val)==seen_positions.end()){
            new_positions.insert(val);  // insert is a no-op if val already exists
        }
    }
    return {new_positions.begin(), new_positions.end()};
}

int main(int argc, char** argv){

    if(argc!=6){
        std::cout<<"usage: ./gen_patterns <input_file> <length> <number> <patterns file> <forbidden symbol>"<<std::endl;
        exit(1);
    }

    std::string input_file = std::string(argv[1]);

    char *pend;
    long int len = strtol(argv[2], &pend, 10);
    long int number = strtol(argv[3], &pend, 10);
    std::string output_file = std::string(argv[4]);

    long int fb_sym_int = strtol(argv[5], &pend, 10);
    assert(fb_sym_int<256);
    auto fb_sym = (char)fb_sym_int;

    auto f_size = (long int)std::filesystem::file_size(input_file);

    number = std::min(f_size-len+1, number);
    len = std::min(f_size, len);
    size_t rem = number;

    std::unordered_set<uint64_t> seen_positions;

    std::string buffer(len, 0);
    std::ifstream ifs(input_file, std::ios::binary);
    std::ofstream ofs(output_file, std::ios::out);

    std::string file_name = std::filesystem::path(input_file).filename();
    std::string header = "# number="+std::to_string(number)+" length="+std::to_string(len)+" file="+file_name+" forbidden="+std::to_string((int)fb_sym)+"\n";
    ofs.write(header.data(), (std::streamsize)header.size());

    while(rem>0){

        std::vector<uint64_t> rd_pos = sample_random_positions(seen_positions, f_size-len+1, rem);
        std::sort(rd_pos.begin(), rd_pos.end());

        for(auto const& pos : rd_pos){
            ifs.seekg((long long)pos);
            ifs.read(buffer.data(), len);
            if(buffer.find(fb_sym)==std::string::npos){
                ofs.write(buffer.data(), len);
                rem--;
            }
            seen_positions.insert(pos);
        }
    }

    ifs.close();
    ofs.close();
    std::cout<<"We extracted "<<number<<" random patterns of length "<<len<<" from file "<<file_name<<std::endl;
}
